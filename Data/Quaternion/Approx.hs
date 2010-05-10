-- |Provide approximate comparison functions within a factor of
-- machine epsilon for Doubles, and an extensible type-class through
-- which this can be implemented for other types.
module Data.Quaternion.Approx (Approx (..), sgn) where

-- |Implement methods for (~==) approx-eq and (~<) approx-lt at the
-- minimum.
class Approx a where
  -- approx-eq
  (~==) :: a -> a -> Bool
  -- approx-less-than
  (~<)  :: a -> a -> Bool
  -- approx-greater-than
  (~>)  :: a -> a -> Bool
  (~>) = flip (~<)
  -- approx-not-eq
  (~/=) :: a -> a -> Bool
  x ~/= y = not (y ~== x)
  -- approx-less-than-or-eq
  (~<=) :: a -> a -> Bool
  x ~<= y = not (y ~< x)
  -- approx-greater-than-or-eq
  (~>=) :: a -> a -> Bool
  (~>=) = flip (~<=)
  -- clamping
  clamp :: (a, a) -> a -> a
  clamp (a,b) x = if x ~< a then a else if b ~< x then b else x
  -- approx-min
  amin :: a -> a -> a
  amin a b = if a ~< b then a else b
  -- approx-max
  amax :: a -> a -> a
  amax a b = if a ~< b then b else a
  -- approx-compare
  acompare :: a -> a -> Ordering
  acompare a b = if a ~< b then LT else if a ~== b then EQ else GT

-- |approx-sgn (-1, 0, or 1 depending on sign of number)
sgn :: (Approx a, Num a) => a -> Int
sgn a = if a ~== 0 then 0 else if a ~< 0 then -1 else 1

-- compute machine epsilon within a factor of 2
doubleEpsilon :: Double
doubleEpsilon = (last . takeWhile (/= 1) . map (+1) . iterate (/2) $ 1) - 1

instance Approx Double where
  a ~== b = abs (a - b) < doubleEpsilon
  a ~< b = a - b < -doubleEpsilon

instance (Approx a, Approx b, Approx c) => Approx (a, b, c) where
  (a1,b1,c1) ~== (a2,b2,c2) = a1 ~== a2 && b1 ~== b2 && c1 ~== c2
  (a1,b1,c1) ~< (a2,b2,c2) = a1 ~< a2 || (a1 ~== a2 && 
                                (b1 ~< b2 || (b1 ~== b2 && c1 ~< c2)))
  clamp ((a1, b1, c1), (a2, b2, c2)) (a3, b3, c3)
    = (clamp (a1, a2) a3, clamp (b1, b2) b3, clamp (c1, c2) c3)
  (a1, b1, c1) `amin` (a2, b2, c2) = (amin a1 a2, amin b1 b2, amin c1 c2)
  (a1, b1, c1) `amax` (a2, b2, c2) = (amax a1 a2, amax b1 b2, amax c1 c2)
