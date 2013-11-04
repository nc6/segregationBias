{-# LANGUAGE GeneralizedNewtypeDeriving,LambdaCase #-}

module Main where
  import Bio.VariantCall.SegregationBias

  import Control.Monad

  import Data.Attoparsec.ByteString.Char8
  import qualified Data.ByteString.Char8 as B
  import qualified Data.Vector.Unboxed as V

  import System.Console.GetOpt
  import System.Environment (getArgs)
  import System.IO

  defaultOptions :: Options
  defaultOptions = Options {
      optVariantSitesAlpha = 1
    , optVariantSitesBeta = 5
    , optTrueVariantRate = 0.95
    , optFalseVariantRateAlpha = 1
    , optFalseVariantRateBeta = 5
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
      , Option [] ["V"] (ReqArg (\n o -> o { optTrueVariantRate = read n }) "P(V)")
          "Prior likelihood on a variant being present at a given site." 
    ]

  usage :: String
  usage = usageInfo header options
    where header = "Calculate segregation bias for a set of samples.\n" ++
                    "Usage: segbias [Option...] template"

  parseSample :: Parser Sample
  parseSample = do
      a <- decimal
      char '\t'
      b <- decimal
      return $ Sample (a, b)
    
  parseSite :: Parser Site
  parseSite = do
    skipWhile (not . (== '\t'))
    char '\t'
    chr <- decimal
    char '\t'
    id' <- decimal
    char '\t'
    samples <- many1 parseSample
    return $ Site (chr) (id') (V.fromList samples)

  processLine :: B.ByteString -> Either String String
  processLine bs = liftM show $ parseOnly parseSite bs

  main :: IO ()
  main = do
    args <- getArgs
    case (getOpt Permute options args) of
      (o,[f],[]) -> processFile f (foldl (flip id) defaultOptions o)
      (_,_,errs) -> putStrLn (concat errs ++ "\n" ++ usage)
  
  processFile :: String -> Options -> IO ()
  processFile fn opts = do
      file <- openFile fn ReadMode
      line <- B.hGetLine file
      let res = processLine line
      case res of
        Right a -> putStrLn a
        Left b -> putStrLn $ "Error: " ++ b