function [iiMax, nP, ZTolerance, YPattern, muFactor] = SH_11_general()
% ------------------------------------------------------------------------------
% 2023-01-21 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% General pattern for symmetric CHM of order N = 11
% with d = 0 and max. possible #L = 3081.
% ------------------------------------------------------------------------------
% Convergence is quite slow and works only for particular initial conditions.
% It should reproduce all other cases.
% ------------------------------------------------------------------------------
% It contains solutions for:
%
% * non-BH with d = 0, #L = 3081 (max. possible)
% ..............................................................................
% * non-BH with d = 0, #L = 1561
% ..............................................................................
% * non-BH with d = 0, #L = 1457
% ..............................................................................
% * non-BH with d = 0, #L = 751
%   ORIGINAL PATTERN:
%   Y =
%   1 1 1 1 1 1 1 1 1 1 1
%   1 a b c d e f g h j k
%   1 b m n o p q r s f t
%   1 c n u v w s x y h z
%   1 d o v . A o B v d A
%   1 e p w A C t D z k E
%   1 f q s o t m r n b p
%   1 g r x B D r . x g D
%   1 h s y v z n x u c w
%   1 j f h d k b g c a e
%   1 k t z A E p D w e C
%
%   with:
%
%   Y(5,5)/A = o
%   Y(5,5)/o = A
%   Y(5,5)/B = B
%   Y(5,5)/v = d
%   Y(5,5)/d = v
%   Y(5,5)/A = o
%   Y(5,5)/A = o
%   Y(5,5)/o = A
%   Y(5,5)/v = d
% ..............................................................................
% * non-BH with d = 0, #L = 425
%   ORIGINAL PATTERN:
%   Y =
%   1 1 1 1 1 1 1 1 1 1 1
%   1 a b c d e f b g b h
%   1 b . j k m n o p b q
%   1 c j . r s t u v w x
%   1 d k r . y z A B C D
%   1 e m s y . E F G H J
%   1 f n t z E . K L M N
%   1 b o u A F K . P b Q
%   1 g p v B G L P . R S
%   1 b b w C H M b R a T
%   1 h q x D J N Q S T .
%
%   with:
%
%   a/b = b
%   a/c = C
%   a/d = M
%   a/e = T
%   a/f = w
%   a/b = b
%   a/g = H
%   a/b = b
%   a/h = R
%   a/b = b
%   a/b = b
%   a/c = C
%   a/w = f
%   a/d = M
%   a/e = T
%   a/H = g
%   a/f = w
%   a/b = b
%   a/b = b
%   a/g = H
%   a/R = h
%   a/b = b
%   a/b = b
%   a/R = h
%   Y(6,6)/G = H
%   Y(6,6)/H = G
%   Y(6,6)/G = H
%   Y(9,9)/R = S
%   Y(9,9)/S = R
%   Y(9,9)/R = S
% ..............................................................................
% * non-BH with d = 0 and #L = 323
%   ORIGINAL PATTERN:
%   Y =
%       1   1   1   1   1   1   1   1   1   1   1  ;
%       1   a   a'  b   b   c'  c'  b'  b'  c   c  ;
%       1   a'  a   c'  c'  b   b   c   c   b'  b' ;
%       1   b   c'  d   h   g   j   e   f   f'  e' ;
%       1   b   c'  h   d   j   g   f   e   e'  f' ;
%       1   c'  b   g   j   d   h   e'  f'  f   e  ;
%       1   c'  b   j   g   h   d   f'  e'  e   f  ;
%       1   b'  c   e   f   e'  f'  j'  g'  d'  h' ;
%       1   b'  c   f   e   f'  e'  g'  j'  h'  d' ;
%       1   c   b'  f'  e'  f   e   d'  h'  j'  g' ;
%       1   c   b'  e'  f'  e   f   h'  d'  g'  j' ;
% ..............................................................................
% * non-BH with d = 0, #L = 191
%   ORIGINAL PATTERN:
%   Y =
%   1 1 1 1 1   1 1 1   1 1   1
%   1 . a b c   d e f   g h   j
%   1 a . k m   n o p   q r   e
%   1 b k . f   s t u   v w   x
%   1 c m f .   y z A   B C   1/k
%   1 d n s y   . D j   E F   G
%   1 e o t z   D . H   J b   K
%   1 f p u A   j H .   L 1/D M
%   1 g q v B   E J L   . N   P
%   1 h r w C   F b 1/D N .   Q
%   1 j e x 1/k G K M   P Q   .
%
%   with (redundancies):
%
%   a/c = k
%   b/j = k
%   b/u = D
%   d/h = D
%   e/f = D
%   e/k = z
%   e/t = 1/k
%   f/k = M
%   f/p = 1/k
%   g/h = e
%   g/j = a
%   g/b = c
%   g/f = d
%   j/D = Q
%   j/K = 1/D
%   Y(3,3)/k = C
%   m/n = f
%   m/o = y
%   m/p = Y(5,5)
%   m/r = A
%   m/e = C
%   m/k = B
%   n/k = Y(5,5)
%   n/C = D
%   n/D = C
%   o/k = A
%   o/D = Y(8,8)
%   p/k = f
%   p/t = 1/D
%   q/r = o
%   q/k = m
%   q/t = C
%   q/n = p
%   q/Y(6,6) = Y(8,8)
%   r/k = y
%   r/Y(6,6) = 1/D
%   s/D = b
%   t/D = p
%   v/w = t
%   v/s = u
%   v/G = p
%   v/b = b
%   w/G = 1/D
%   x/f = G
%   x/s = Q
%   x/j = b
%   x/e = w
%   x/t = Y(11,11)
%   x/K = u
%   x/v = 1/k
%   m/f = n
%   m/Y(5,5) = p
%   y/D = Y(10,10)
%   z/D = M
%   A/Y(7,7) = 1/D
%   B/C = z
%   B/n = M
%   B/y = A
%   B/Y(7,7) = Y(10,10)
%   d/D = h
%   n/D = C
%   s/D = b
%   y/D = Y(10,10)
%   Y(6,6)/D = r
%   j/b = 1/k
%   E/F = D
%   E/K = r
%   F/G = C
%   F/D = N
%   F/j = r
%   G/D = w
%   e/t = 1/k
%   e/D = f
%   o/D = Y(8,8)
%   t/D = p
%   z/D = M
%   H/J = 1/D
%   H/b = A
%   H/K = Y(8,8)
%   H/f = e
%   H/A = b
%   H/Y(8,8) = K
%   H/1/D = J
%   f/1/D = e
%   u/1/D = b
%   u/b = 1/D
%   j/1/D = K
%   j/b = 1/k
%   H/Y(8,8) = K
%   g/h = e
%   N/r = Q
%   N/w = C
%   N/b = Y(10,10)
%   P/Q = K
%   P/e = Y(11,11)
%   P/x = 1/k
%   P/G = M
% ..............................................................................
% * non-BH with d = 0, #L = 63
%   ORIGINAL PATTERN
%   Y =
%   1 1 1 1 1 1 1   1 1   1   1
%   1 a b c d e f   g h   j   k
%   1 b a g d e k   c j   h   f
%   1 c g m n f k   c f   k   f
%   1 d d n . o n   n d   d   n
%   1 e e f o . k   f p   p   k
%   1 f k k n k q   f r   1/g r
%   1 g c c n f f   m k   f   k
%   1 h j f d p r   k s   t   1/g
%   1 j h k d p 1/g f t   s   r
%   1 k f f n k r   k 1/g r   q
%
%   with:
%
%   a/b = k
%   a/c = h
%   a/d = g
%   a/e = e
%   a/f = j
%   b/g = j
%   b/j = g
%   c/f = g
%   c/g = f
%   d/e = p
%   d/f = h
%   d/g = s
%   d/j = k
%   d/b = r
%   d/k = j
%   d/c = t
%   e/g = p
%   e/p = g
%   f/g = k
%   h/g = t
%   j/h = g
%   j/g = h
%   k/g = r
%   b/g = j
%   k/r = g
%   m/n = g
%   m/f = f
%   m/k = c
%   n/f = k
%   n/c = r
%   n/q = g
%   Y(5,5)/o = o
%   Y(5,5)/n = d
%   Y(5,5)/m = s
%   o/e = k
%   o/f = p
%   q/f = r
%   q/k = k
%   s/t = f
%   s/j = r
%   s/h = k
% ------------------------------------------------------------------------------

    iiMax = 40000;
    nP = 45;
    ZTolerance = 6e-13;
    YPattern = @pattern;
    muFactor = 0.0007;
end

function Y = pattern(p)
    x1 = -(1+p(1)+p(2)+p(3)+p(4)+p(5)+p(6)+p(7)+p(8)+p(9));
    x2 = -(1+p(1)+p(10)+p(11)+p(12)+p(13)+p(14)+p(15)+p(16)+p(17));
    x3 = -(1+p(2)+p(10)+p(18)+p(19)+p(20)+p(21)+p(22)+p(23)+p(24));
    x4 = -(1+p(3)+p(11)+p(18)+p(25)+p(26)+p(27)+p(28)+p(29)+p(30));
    x5 = -(1+p(4)+p(12)+p(19)+p(25)+p(31)+p(32)+p(33)+p(34)+p(35));
    x6 = -(1+p(5)+p(13)+p(20)+p(26)+p(31)+p(36)+p(37)+p(38)+p(39));
    x7 = -(1+p(6)+p(14)+p(21)+p(27)+p(32)+p(36)+p(40)+p(41)+p(42));
    x8 = -(1+p(7)+p(15)+p(22)+p(28)+p(33)+p(37)+p(40)+p(43)+p(44));
    x9 = -(1+p(8)+p(16)+p(23)+p(29)+p(34)+p(38)+p(41)+p(43)+p(45));
    x10= -(1+p(9)+p(17)+p(24)+p(30)+p(35)+p(39)+p(42)+p(44)+p(45));

    Y = [
        1   1       1       1       1       1       1       1       1       1       1      ;
        1   x1      p(1)    p(2)    p(3)    p(4)    p(5)    p(6)    p(7)    p(8)    p(9)   ;
        1   p(1)    x2      p(10)   p(11)   p(12)   p(13)   p(14)   p(15)   p(16)   p(17)  ;
        1   p(2)    p(10)   x3      p(18)   p(19)   p(20)   p(21)   p(22)   p(23)   p(24)  ;
        1   p(3)    p(11)   p(18)   x4      p(25)   p(26)   p(27)   p(28)   p(29)   p(30)  ;
        1   p(4)    p(12)   p(19)   p(25)   x5      p(31)   p(32)   p(33)   p(34)   p(35)  ;
        1   p(5)    p(13)   p(20)   p(26)   p(31)   x6      p(36)   p(37)   p(38)   p(39)  ;
        1   p(6)    p(14)   p(21)   p(27)   p(32)   p(36)   x7      p(40)   p(41)   p(42)  ;
        1   p(7)    p(15)   p(22)   p(28)   p(33)   p(37)   p(40)   x8      p(43)   p(44)  ;
        1   p(8)    p(16)   p(23)   p(29)   p(34)   p(38)   p(41)   p(43)   x9      p(45)  ;
        1   p(9)    p(17)   p(24)   p(30)   p(35)   p(39)   p(42)   p(44)   p(45)   x10    ;
    ];
end
