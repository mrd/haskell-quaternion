-- |Quaternion (a.k.a. 4-dimensional complex numbers) operations, with
-- particular emphasis on the representation of rotations in
-- 3-dimensional space.
--
-- Typical quaternion usage in computer graphics is to construct or
-- compose one, such as a 45-degree rotation around the X-axis:
--   @let q = angleAxis 45 (1, 0, 0)@
-- and then either directly multiply it by a vector like so:
--   @toVec (reciprocal q `mul` fromVec v `mul` q)@
-- or perhaps convert it to a matrix for use with OpenGL:
--   @makeMatrix (4, 4) (rowMajorElems q)@

module Data.Quaternion
  ( Quat, zero, real, imag, mul
  , rowMajorElems, eulerAngles, angleAxis
  , normalize, reciprocal, slerp, lookAt
  , fromVec, toVec, dot, plus, scale ) 
where

data Quat a = Q { real :: a, imag :: (a, a, a) }
type Vector a = (a, a, a)

-- |The Zero quaternion
zero :: Num a => Quat a
zero = Q { real = 1, imag = (0, 0, 0) }

-- |Quaternion dot product
dot :: (Num a) => Quat a -> Quat a -> a
q `dot` r = real q * real r + dotProduct (imag q) (imag r)

-- |Quaternion addition
plus :: (Num t) => Quat t -> Quat t -> Quat t
q `plus` r = Q { real = real q + real r
               , imag = (qx + rx, qy + ry, qz + rz) }
  where (qx, qy, qz) = imag q
        (rx, ry, rz) = imag r

-- |Scale quaternion by constant factor
scale :: (Num t) => t -> Quat t -> Quat t
n `scale` q = Q { real = real q * n, imag = (qx * n, qy * n, qz * n) }
  where (qx, qy, qz) = imag q

-- |Quaternion multiplication
mul :: (Num t) => Quat t -> Quat t -> Quat t
q1 `mul` q2 = Q { real = s1 * s2 - dotProduct (imag q1) (imag q2)
                , imag = (x4, y4, z4) }
  where
    (s1, s2)     = (real q1, real q2)
    (x1, y1, z1) = imag q1
    (x2, y2, z2) = imag q2
    (x3, y3, z3) = crossProduct (imag q1) (imag q2)
    x4           = x3 + s1 * x2 + s2 * x1
    y4           = y3 + s1 * y2 + s2 * y1
    z4           = z3 + s1 * z2 + s2 * z1

-- |Conversion to 4x4 homogeneous matrix elements in row-major order
rowMajorElems :: (Num t) => Quat t -> [t]
rowMajorElems q = es
  where
    s         = real q
    (a, b, c) = imag q
    es = [ 1 - 2*b*b - 2*c*c, 2*a*b + 2*s*c     , 2*a*c - 2*s*b,     0
         , 2*a*b - 2*s*c    , 1 - 2*a*a - 2*c*c , 2*b*c + 2*s*a,     0
         , 2*a*c + 2*s*b    , 2*b*c - 2*s*a     , 1 - 2*a*a - 2*b*b, 0
         , 0                , 0                 , 0                , 1 ]

-- |Quaternion normalization to unit length
normalize :: (Floating t) => Quat t -> Quat t
normalize q = Q { real = s / mag, imag = (a / mag, b / mag, c / mag) }
  where
    s         = real q
    (a, b, c) = imag q
    mag = sqrt (s*s + a*a + b*b + c*c)

-- |Quaternion reciprocal
reciprocal :: (Fractional a) => Quat a -> Quat a
reciprocal q = Q { real = real q / mag2, imag = (-x / mag2, -y / mag2, -z / mag2) }
  where (x,y,z) = imag q
        s = real q
        mag2 = (s*s+x*x+y*y+z*z)

-- |Given an angle and an axis, produce a quaternion describing that
-- rotation.
angleAxis :: (Floating a) => 
                 a        -- ^ Angle in degrees.
              -> Vector a -- ^ Axis of rotation
              -> Quat a
angleAxis a (x, y, z) = Q { real = c, imag = (s * x', s * y', s * z') }
  where
    s = sin (pi * a / 360)
    c = cos (pi * a / 360)
    (x', y', z') = (x / mag, y / mag, z / mag)
    mag = sqrt (x*x+y*y+z*z)

-- |Compute a quaternion that rotates from the current direction to
-- the desired direction, with a given notion of "up".
lookAt :: (Floating t, Ord t) => 
             Vector t -- ^ Current direction of object
          -> Vector t -- ^ Desired direction of object
          -> Vector t -- ^ Vector considered to be "up"
          -> Quat t
lookAt cur tgt up 
  | nx == 0 && ny == 0 && nz == 0 = 
    -- the normal magnitude is 0 so the vectors are parallel
    if sgn cx == sgn tx && sgn cy == sgn ty && sgn cz == sgn tz then
    -- sgn of cur is the same as tgt means they point in the same dir
      Q { real = 1, imag = (0, 0, 0) }
    -- doing a 180-degree rotation, be careful of axis
    else if cx == 0 then
      angleAxis 180 (crossProduct cur (1, 0, 0))
    else
      angleAxis 180 (crossProduct cur (0, 0, 1))
  -- if projected up-vector is meaningless then just forget trying to re-orient
  | up_p == (0, 0, 0) = q 
  | otherwise          = q' `mul` q
  where   
    (cx, cy, cz) = cur
    (tx, ty, tz) = tgt
    (ux, uy, uz) = up
    (nx, ny, nz) = crossProduct cur tgt
    angle        = acos' ( dotProduct cur tgt /
                           sqrt ((cx*cx+cy*cy+cz*cz) * (tx*tx+ty*ty+tz*tz)))
    q    = angleAxis (r2d angle) (nx, ny, nz)
    upq  = Q { real = 0, imag = up }
    upq' = q `mul` upq `mul` reciprocal q
    -- up' is the up-vector rotated according to q
    up'  = imag upq'
    -- project up and up' onto the plane with normal tgt
    up_p = project up tgt
    up_p' = project up' tgt
    mags  = dotProduct up_p up_p * dotProduct up_p' up_p'
    -- find angle between them and rotate around normal
    xprod = crossProduct up_p' up_p
    beta  = acos' ( dotProduct up_p up_p' / sqrt mags )
    -- first obtain magnitude of rotation, then determine by
    -- cross-product if the rotation is in the correct direction.
    beta' = if sameDir xprod tgt then beta else -beta
    q'    = angleAxis (r2d beta') tgt

-- |Convert the quaternion into an equivalent triple of Euler angles.
eulerAngles :: (RealFloat t) => Quat t -> Vector t
eulerAngles q = ( r2d (atan2 (2*s*a + 2*b*c) (1 - 2*(a*a + b*b)))
                , r2d (asin (2*(s*b - c*a)))
                , r2d (atan2 (2*s*c+2*a*b) (1 - 2*(b*b + c*c))) )
  where
    s         = real q
    (a, b, c) = imag q

-- |Use the SLERP method of interpolating a partial rotation between
-- two quaternions.
slerp :: (Floating a, Ord a) => 
           Quat a -- ^ Starting quaternion
        -> Quat a -- ^ Ending quaternion
        -> a      -- ^ Interpolation factor
        -> Quat a
slerp q0 q1 t 
  -- if inputs are not co-linear, use SLERP  
  | d < 1    = (cos theta `scale` nq0) `plus` (sin theta `scale` nq2)
  -- inputs are too close to being co-linear, so just interpolate  
  | otherwise = normalize (nq0 `plus` (t `scale` (nq1 `plus` ((-1) `scale` nq0))))
  where
    nq0 = normalize q0
    nq1 = normalize q1
    d = nq0 `dot` nq1
    d' = min 1 (max (-1) d) -- clamp [-1, 1]
    theta = t * acos d
    q2 = nq1 `plus` ((-d) `scale` nq0)
    nq2 = normalize q2

-- |Produce a zero-rotation quaternion from a vector
fromVec :: (Num a) => Vector a -> Quat a
fromVec v = Q { real = 0, imag = v }

-- |Discard the real part and return a vector
toVec :: Quat a -> Vector a
toVec = imag

--------------------------------------------------
-- Utility

normv (x,y,z) = (check0 x/mag,check0 y/mag,check0 z/mag)
  where mag = sqrt (x*x+y*y+z*z)

sameDir (x1,y1,z1) (x2,y2,z2) = 
  sgn x1 == sgn x2 &&
  sgn y1 == sgn y2 &&
  sgn z1 == sgn z2

check0 x = if x == 0 then 0 else x
check1 x = if x == 1 then 1 else x

-- prevent minor floating point errors from screwing up acos
acos' x = acos (check0 (check1 x))

crossProduct (a0, a1, a2) (b0, b1, b2) =
  (a1*b2 - a2*b1, a2*b0 - a0*b2, a0*b1 - a1*b0)

dotProduct (a0, a1, a2) (b0, b1, b2) = a0*b0 + a1*b1 + a2*b2

r2d a = 180 * a / pi

project (vx, vy, vz) (nx, ny, nz) = (ux, uy, uz)
  where
    nmag          = sqrt (nx*nx + ny*ny + nz*nz)
    (nx',ny',nz') = (nx / nmag, ny / nmag, nz / nmag)
    dp            = vx * nx' + vy * ny' + vz * nz'
    ux            = vx - nx' * dp
    uy            = vy - ny' * dp
    uz            = vz - nz' * dp

sgn x | x < 0 = -1 | x > 0 = 1 | x == 0 = 0
