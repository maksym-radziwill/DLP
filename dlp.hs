import GHC.Integer.GMP.Internals
import Data.List
import System.Environment

---- Quick implementation of Pollard Rho
---- One could get away with something much more primitive

f :: Integer -> (Integer -> Integer -> Integer)
f k = g k

g :: Integer -> Integer -> Integer -> Integer      
g k x n = mod (x^k + 1) n

pollard :: Integer -> Integer -> Integer -> Integer -> (Integer -> Integer -> Integer) -> Integer
pollard x y d n g | d == 1 =
                    let x0 = g x n
                        y0 = g (g y n) n
                    in pollard x0 y0 (gcd (abs (x0 - y0)) n) n g
                  | d == n = 0
                  | otherwise = d

factor :: Integer -> [Integer] -> [Integer]                                
factor n primelist = let k x = (mod n x) == 0
                         divlist = filter k primelist
                         prod = product divlist
                     in filter (/= 1) (factor' (div n prod) (1:divlist) primelist)
    where
      factor' n fac primelist | (nextPrimeInteger (n-1)) == n = n : fac
                              | n == 1 = fac
                              | otherwise =
                                  let pol = [pollard y x 1 n (f k) | k <- [2..] , x <- [2..5] , y <- [2..5] ]
                                      d = head (filter (/= 0) pol)
                                  in factor' (div n d) ((factor' d [1] primelist) ++ fac) primelist


factorization :: (Eq a) => [a] -> [(a, Integer)]
factorization fac | null fac = []
                  | otherwise =
                      let h = head fac
                      in ( h, freq h fac ) : (factorization (filter (/= h) fac)) 

factors :: Integer -> [Integer] -> [(Integer, Integer)]                         
factors n primelist = factorization (factor n primelist)
                         
freq :: (Eq a) => a -> [a] -> Integer
freq x s | s == [] = 0
         | otherwise = if x == head s
                       then 1 + (freq x (tail s))
                       else (freq x (tail s))

genprimes :: Integer -> [Integer]
genprimes n = genprimes' 1 n
    where genprimes' k n | k > n = []
                         | otherwise = let p = nextPrimeInteger k
                                       in p : (genprimes' p n)

---- Chinese remainder theorem borrowed from Rosetta code because I was too lazy to write it
              
egcd :: Integral a => a -> a -> (a, a)
egcd _ 0 = (1, 0)
egcd a b = (t, s - q * t)
  where
    (s, t) = egcd b r
    (q, r) = a `quotRem` b
 
modInv :: Integral a => a -> a -> Maybe a
modInv a b = case egcd a b of
  (x, y) | a * x + b * y == 1 -> Just x
         | otherwise          -> Nothing
 
chineseRemainder :: Integral a => [a] -> [a] -> Maybe a
chineseRemainder residues modulii = do
  inverses <- sequence $ zipWith modInv crtModulii modulii
  return . (`mod` modPI) . sum $
    zipWith (*) crtModulii (zipWith (*) residues inverses)
  where
    modPI = product modulii
    crtModulii = map (modPI `div`) modulii

unJust :: Maybe a -> a
unJust (Just x) = x

crt :: [(Integer, Integer)] -> Integer
crt xs = unJust (chineseRemainder [ fst x | x <- xs ] [ snd x | x <- xs ] )

---- Primitive roots

checkPrimRoot :: Integer -> Integer -> Bool               
checkPrimRoot g p = all (/= 1) [ fastPow (g `mod` p) (div (p-1) q) p | q <- (map fst (factors (p-1) (genprimes 4000))) ]
         
----- DLP for p such that p - 1 has small square-free prime factors
         
dlp :: Integer -> Integer -> Integer -> Integer
dlp p g h = crt [ (bruteDlp p (fastPow g (div (p-1) (q^k)) p) (fastPow h (div (p-1) (q^k)) p), q^k) | (q,k) <- (factors (p-1) (genprimes 4000)) ]

bruteDlp :: Integer -> Integer -> Integer -> Integer
bruteDlp p g h = bruteDlp' p g h 1
    where bruteDlp' p0 g0 h0 x0 = let v = fastPow g0 (mod x0 (p0 - 1)) p0
                                  in if ((mod v p0) == (mod h0 p0))
                                     then x0
                                     else bruteDlp' p0 g0 h0 (x0 + 1)

---- Fast exponentiation

fastPow :: Integer -> Integer -> Integer -> Integer
fastPow 1 k p = 1                                          
fastPow g 0 p = 1                                          
fastPow g k p = let a = (g^(mod k 2)) `mod` p
                    v = mod (g*g) p 
                in mod (a * (fastPow v (div (k - (mod k 2)) 2) p)) p                                      

---- Main loop
        
toInt :: String -> Integer
toInt x = read x :: Integer

main :: IO ()
main = do
  args <- getArgs
  if ((length args) /= 3)
  then print "Usage : prog p g h"
  else do
    [p,g,h] <- getArgs
    if (nextPrimeInteger ((toInt p) - 1) /= (toInt p))
    then print "Not prime!"
    else do
      if (checkPrimRoot (toInt g) (toInt p))
      then print (dlp (toInt p) (toInt g) (toInt h))
      else print "Not a primitive root"

    
