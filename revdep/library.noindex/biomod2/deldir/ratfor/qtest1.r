subroutine qtest1(h,i,j,k,x,y,ntot,eps,shdswp,nerror)

# The Lee-Schacter test for the LOP (all points are real,
# i.e. non-ideal).  If the LOP is ***not*** satisfied (i.e. if
# vertex j is inside the circumcircle of vertices h, i, and k) then the
# diagonals should be swapped, i.e. shdswp ("should-swap") is true.
# Called by qtest.

implicit double precision(a-h,o-z)
dimension x(-3:ntot), y(-3:ntot)
integer h
logical shdswp

# The vertices of the quadrilateral are labelled
# h, i, j, k in the anticlockwise direction, h
# being the point of central interest.

# Make sure the quadrilateral is convex, so that
# it makes sense to swap the diagonal.
# call acchk(i,j,k,shdswp,x,y,ntot,eps)
# if(!shdswp) return
#
# 23 July 2011:
# The foregoing test is a load of dingos' kidneys.  (1) It is
# unnecessary, and (2) it is wrong!  (1) If the LOP is not satisfied
# (the only circumstance under which there should be a swap) then the
# quadrilateral ***must*** be convex, and so swapping can sensibly
# take place.  (2) The vertices i, j, k in will ***always*** be in
# anticlockwise order, since the vertices h, i, j, k of the quadrilateral
# are in such order and i is connected to k, whence j can't be inside
# the triangle ihk.  So the test does nothing.  But then it didn't need
# to do anything.

# Get the coordinates of vertices h and j.
xh = x(h)
yh = y(h)
xj = x(j)
yj = y(j)

# Find the centre of the circumcircle of vertices h, i, k.
call circen(h,i,k,x0,y0,x,y,ntot,eps,shdswp,nerror)
if(nerror>0) return
if(shdswp) return # The points h, i, and k are colinear, so
                  # the circumcircle has `infinite radius', so
                  # (xj,yj) is definitely inside.

# Check whether (xj,yj) is inside the circle of centre
# (x0,y0) and radius r = dist[(x0,y0),(xh,yh)]

a  = x0-xh
b  = y0-yh
r2 = a*a+b*b
a  = x0-xj
b  = y0-yj
ch = a*a + b*b
if(ch<r2) shdswp = .true.
else shdswp = .false.

return
end
