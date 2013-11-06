{-# LANGUAGE GeneralizedNewtypeDeriving,LambdaCase #-}

module Main where
  import Bio.VariantCall.SegregationBias

  import Control.Applicative (optional)
  import Control.Monad
  import Control.Monad.IO.Class

  import Data.Attoparsec.ByteString.Char8
  import qualified Data.ByteString.Char8 as B
  import qualified Data.Vector.Unboxed as V
  import Data.Conduit
  import qualified Data.Conduit.Binary as CB
  import qualified Data.Conduit.List as CL

  --import Debug.Trace (traceShow)

  import System.Console.GetOpt
  import System.Environment (getArgs)

  defaultOptions :: Options
  defaultOptions = Options {
      optVariantSitesAlpha = 1
    , optVariantSitesBeta = 4
    , optTrueVariantRate = 0.95
    , optFalseVariantRateAlpha = 1
    , optFalseVariantRateBeta = 4
    , optVariantFreq = 0.5
  }

  options :: [OptDescr (Options -> Options)]
  options =
    [
        Option [] ["M_alpha"] (ReqArg (\n o -> o { optVariantSitesAlpha = read n }) "ALPHA")
          "Alpha parameter for beta-binomial governing P(M|V)."
      , Option [] ["M_beta"] (ReqArg (\n o -> o { optVariantSitesBeta = read n }) "BETA")
          "Beta parameter for beta-binomial governing P(M|V)."
      , Option [] ["p_T"] (ReqArg (\n o -> o { optTrueVariantRate = read n }) "p_T")
          "Probability of seeing a variant call having observed a variant at a given site."
      , Option [] ["P_F_alpha"] (ReqArg (\n o -> o { optFalseVariantRateAlpha = read n }) "p_F_alpha")
          "Alpha parameter for the Beta distributed probability of seeing a variant call when no variant is present."
      , Option [] ["P_F_beta"] (ReqArg (\n o -> o { optFalseVariantRateAlpha = read n }) "p_F_beta")
          "Beta parameter for the Beta distributed probability of seeing a variant call when no variant is present."
      , Option [] ["V"] (ReqArg (\n o -> o { optVariantFreq = read n }) "P(V)")
          "Prior likelihood on a variant being present at a given site." 
    ]

  usage :: String
  usage = usageInfo header options
    where header = "Calculate segregation bias for a set of samples.\n" ++
                    "Usage: segbias [Option...] template"

  parseSample :: Parser Sample
  parseSample = do
      d <- decimal
      _ <- char '\t'
      o <- decimal
      _ <- optional $ char '\t'
      return $ Sample (o, o+d)
    
  parseSite :: Parser Site
  parseSite = do
    skipWhile (not . (== '\t'))
    _ <- char '\t'
    chr <- decimal
    _ <- char '\t'
    id' <- decimal
    _ <- char '\t'
    samples <- many1 parseSample
    return $ Site (chr) (id') (V.fromList samples)

  processLine :: Options -> B.ByteString -> Either String Double
  processLine opts bs = liftM go $ parseOnly parseSite bs
    where go (Site _ _ samples) = probVariantGivenObs opts samples

  main :: IO ()
  main = do
    args <- getArgs
    case (getOpt Permute options args) of
      (o,[f],[]) -> processFile f (foldl (flip id) defaultOptions o)
      (_,_,errs) -> putStrLn (concat errs ++ "\n" ++ usage)
  
  processFile :: String -> Options -> IO ()
  processFile fn opts = let
      source = CB.sourceFile fn $= CB.lines
      process = CL.map $ processLine opts
      sink = CL.mapM_ $ liftIO . showResult
      showResult = \case
        Right a -> putStrLn $ show a
        Left b -> putStrLn $ "Error: " ++ b
    in runResourceT $ source $= process $$ sink