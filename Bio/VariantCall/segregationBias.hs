{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module Bio.VariantCall.SegregationBias (
      probVariantGivenObs
    , Options(..)
    , Sample(..)
    , Site(..)
  ) where

  import Data.Vector.Unboxed (Vector)
  import qualified Data.Vector.Unboxed as V
  import qualified Data.Vector.Generic.Base
  import qualified Data.Vector.Generic.Mutable

  --import Debug.Trace (traceShow, trace)

  import qualified Numeric.Integration.TanhSinh as I (absolute, parTrap, result)
  import Numeric.SpecFunctions (choose, logBeta)

  import qualified Statistics.Distribution as D
  import Statistics.Distribution.BetaBinomial (betaBinomialDistr)
  import Statistics.Distribution.Binomial (binomial)

  -- | Sample (Variant, total)
  newtype Sample = Sample (Int, Int) deriving (
      Eq
      , Show
      , Data.Vector.Generic.Base.Vector Vector
      , Data.Vector.Generic.Mutable.MVector V.MVector
      , V.Unbox
    )

  -- | Site chr id samples
  data Site = Site Int Int (Vector Sample) deriving (Eq, Show)

  -- | Parameters for the calculation.
  data Options = Options {
      optVariantSitesAlpha :: Double -- ^ Alpha parameter for beta-binomial governing P(M|V)
    , optVariantSitesBeta :: Double -- ^ Beta parameter
    , optTrueVariantRate :: Double -- ^ p_T
    , optFalseVariantRateAlpha :: Double -- ^ alpha_F
    , optFalseVariantRateBeta :: Double -- ^ beta_F
    , optVariantFreq :: Double -- ^ P(V)
  }

  -- Define some type aliases for notational convenience
  type AlleleCount = Int
  type Probability = Double

  probObsVariantHNR :: Probability -- ^ p_T
          -> AlleleCount -- ^ Total allele count (d)
          -> AlleleCount -- ^ Variant allele count (o)
          -> Probability -- ^ Probability of observing a variant allele given that 
                         --   the individual is HNR
  probObsVariantHNR p_T d o = D.probability binDist o where
    binDist = binomial d p_T

  probObsVariantHet :: Probability -- ^ p_T
          -> Double -- ^ alpha_F
          -> Double -- ^ beta_F
          -> AlleleCount -- ^ Total allele count (d)
          -> AlleleCount -- ^ Variant allele count (o)
          -> Probability -- ^ Probability of observing a variant allele given that the individual is het
  probObsVariantHet _ _ _ 0 0 = 1
  probObsVariantHet _ _ _ 0 _ = 0
  probObsVariantHet p_T a_F b_F d o = result where
    result = (d `choose` o) * int / ((2 ^ d) * (exp $ logBeta a_F b_F))
    o' = fromIntegral o
    d' = fromIntegral d
    int = I.result . I.absolute 1e-6 $ I.parTrap fun 0 1
    fun p_F =    (p_T + p_F)**o'
               * (1 - p_T - p_F)**(d' - o')
               * (p_F)**(a_F -1)
               * (1-p_F)**(b_F -1)

  -- | Calculate the probability of observing a variant allele at a given site given that the indivudual
  -- is homogenous reference.
  probObsVariantHR :: Double -- ^ alpha_F
           -> Double -- ^ beta_F
           -> AlleleCount -- ^ Total allele count (d)
           -> AlleleCount -- ^ Variant allele count (o)
           -> Probability -- ^ Probability of observing a variant allele given that the individual is homogenous ref
  probObsVariantHR a_F b_F d o = D.probability bbDist o where
    bbDist = betaBinomialDistr a_F b_F d

  probObsGivenVariantSite :: Probability -- ^ p_M
            -> Probability -- ^ p_T
            -> Double -- ^ alpha_F
            -> Double -- ^ beta_F
            -> AlleleCount -- ^ Total allele count (d)
            -> AlleleCount -- ^ Variant allele count (o)
            -> Probability
  probObsGivenVariantSite pm p_t a_F b_F d o = let
      het = probObsVariantHet p_t a_F b_F d o
      hnr = probObsVariantHNR p_t d o
    in 2*(1-pm)*het + pm*hnr

  probMGivenV :: Int -- ^ Number of samples (N)
              -> Int -- ^ Number of samples showing variant allele (M)
              -> Double -- ^ alpha
              -> Double -- ^ beta
              -> Probability
  probMGivenV 0 0 _ _ = 1
  probMGivenV 0 _ _ _ = 0
  probMGivenV n m a b = D.probability bbDist m where
    bbDist = betaBinomialDistr a b n

  -- | P(o|V)
  probObsGivenVariant :: Options
                      -> Vector Sample
                      -> Probability
  probObsGivenVariant opts samples = sum foo where
    n = V.length samples
    p_t = optTrueVariantRate opts
    a_F = optFalseVariantRateAlpha opts
    b_F = optFalseVariantRateBeta opts
    a = optVariantSitesAlpha opts
    b = optVariantSitesBeta opts
    foo = map (\m -> (probMGivenV n m a b)*(bar m)) [0 .. n]
    bar m = V.product $ V.map (\(Sample (o,d)) -> baz m d o) samples
    baz m d o = let pm = (fromIntegral m) / (2*(fromIntegral n)) in
                pm * (probObsGivenVariantSite pm p_t a_F b_F d o) + 
                ((1-pm)^2) * (probObsVariantHR a_F b_F d o)

  -- | P(o|!V)
  probObsGivenNonVariant :: Options 
                         -> Vector Sample
                         -> Probability
  probObsGivenNonVariant opts samples = V.product x where
    x = V.map pSite samples
    bbd n = betaBinomialDistr (optFalseVariantRateAlpha opts) (optFalseVariantRateBeta opts) n
    pSite (Sample (o,d)) = D.probability (bbd d) o

  -- | Probability that this is a true variant given the observations.
  probVariantGivenObs :: Options
                      -> Vector Sample
                      -> Probability
  probVariantGivenObs opts samples = 
   (prior * likelihood) / probability where
      prior = optVariantFreq opts
      likelihood = probObsGivenVariant opts samples
      unlikelihood = probObsGivenNonVariant opts samples
      probability = (likelihood * prior) + (1 - prior)*unlikelihood
