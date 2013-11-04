{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module Bio.VariantCall.SegregationBias where

  import Data.Vector.Unboxed (Vector)
  import qualified Data.Vector.Unboxed as V
  import qualified Data.Vector.Generic.Base
  import qualified Data.Vector.Generic.Mutable

  import qualified Statistics.Distribution as D
  import Statistics.Distribution.BetaBinomial (betaBinomialDistr)
  import Statistics.Distribution.Binomial (binomial)

  newtype Sample = Sample (Int, Int) deriving (
    Eq
    , Show
    , Data.Vector.Generic.Base.Vector Vector
    , Data.Vector.Generic.Mutable.MVector V.MVector
    , V.Unbox
    )

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
          -> Probability -- ^ Probability of observing a variant allele given that the individual is probObsVariantHNR
  probObsVariantHNR p_T d o = D.probability binDist o where
    binDist = binomial d p_T

  probObsVariantHet :: Probability -- ^ p_T
          -> Double -- ^ alpha_F
          -> Double -- ^ beta_F
          -> AlleleCount -- ^ Total allele count (d)
          -> AlleleCount -- ^ Variant allele count (o)
          -> Probability -- ^ Probability of observing a variant allele given that the individual is het
  probObsVariantHet = undefined

  -- | Calculate the probability of observing a variant allele at a given site given that the indivudual
  -- is homogenous reference.
  probObsVariantHR :: Double -- ^ alpha_F
           -> Double -- ^ beta_F
           -> AlleleCount -- ^ Total allele count (d)
           -> AlleleCount -- ^ Variant allele count (o)
           -> Probability -- ^ Probability of observing a variant allele given that the individual is homogenous ref
  probObsVariantHR a_F b_F d o = D.probability bbDist o where
    bbDist = betaBinomialDistr a_F b_F d

  probMGivenV :: Int -- ^ Number of samples (N)
              -> Int -- ^ Number of samples showing variant allele (M)
              -> Double -- ^ alpha
              -> Double -- ^ alpha
              -> Probability
  probMGivenV = undefined

  probObsGivenVariant :: Options
                      -> Vector Sample
                      -> Probability
  probObsGivenVariant = undefined

  probObsGivenNonVariant :: Options 
                         -> Vector Sample
                         -> Probability
  probObsGivenNonVariant opts samples = V.product x where
    x = V.map pSite samples
    bbd n = betaBinomialDistr (optFalseVariantRateAlpha opts) (optFalseVariantRateBeta opts) n
    pSite (Sample (d,o)) = D.probability (bbd d) o

  probVariantGivenObs :: Options
                      -> Vector Sample
                      -> Probability
  probVariantGivenObs opts samples = (prior * likelihood) / probability where
    prior = optVariantFreq opts
    likelihood = probObsGivenVariant opts samples
    probability = (likelihood * prior) + (1 - prior)*(probObsGivenNonVariant opts samples)
