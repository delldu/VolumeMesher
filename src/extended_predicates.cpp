#include "extended_predicates.h"
#include <algorithm>

//----------------------
// Derivated predicates
//----------------------

// ----- Class1: the highest dimensional object is a POINT -------

//  Input: point p=(px,py,pz), q=(qx,qy,qz), r=(rx,ry,rz).
// Output: 1 -> if the 3 input points are NOT aligned,
//         0 -> otherwise
bool misAlignment(const double *p, const double *q, const double *r) {

  // Projection on (x,y)-plane
  if (orient2d(p, q, r))
    return 1;

  // Projection on (y,z)-plane
  if (orient2d(p + 1, q + 1, r + 1))
    return 1;

  // Projection on (x,z)-plane
  const double pxz[] = {p[0], p[2]};
  const double qxz[] = {q[0], q[2]};
  const double rxz[] = {r[0], r[2]};
  return (orient2d(pxz, qxz, rxz));
}

// ----- Class2: the highest dimensional object is a SEGMENT -------

//  Input: point p=(px,py,pz); segment v1-v2 with v1=(v1x,v1y,v1z),
//  v2=(v2x,v2y,v2z).
// Output: 1 -> if the point p belong to the interior of the segment v1-v2
// (endpoints excluded),
//         0 -> otherwise.
uint32_t pointInInnerSegment(const double *p, const double *v1,
                             const double *v2) {
  return (misAlignment(p, v1, v2) == 0 &&
          ((v1[0] < v2[0] && v1[0] < p[0] && p[0] < v2[0]) ||
           (v1[0] > v2[0] && v1[0] > p[0] && p[0] > v2[0]) ||
           (v1[1] < v2[1] && v1[1] < p[1] && p[1] < v2[1]) ||
           (v1[1] > v2[1] && v1[1] > p[1] && p[1] > v2[1]) ||
           (v1[2] < v2[2] && v1[2] < p[2] && p[2] < v2[2]) ||
           (v1[2] > v2[2] && v1[2] > p[2] && p[2] > v2[2])));
}

//  Input: point p=(px,py,pz); segment v1-v2 with v1=(v1x,v1y,v1z),
//  v2=(v2x,v2y,v2z).
// Output: 1 -> if the point p belong to the segment v1-v2 (endpoints included),
//         0 -> otherwise.
uint32_t pointInSegment(const double *p, const double *v1, const double *v2) {
  return (same_point(p, v1) || same_point(p, v2) || pointInInnerSegment(p, v1, v2));
}

//  Input: point p, point q, and a segment v1-v2, througt their coordinates:
//         p=(px,py,pz), q=(qx,qy,qz), v1=(v1x,v1y,v1z), v2=(v2x,v2y,v2z).
// Output: 1 -> points p and q lies both on the same side of the straight line
// passing througt v1 and v2.
//         0 -> otherwise.
// Note. Points and segment must be complanar.
uint32_t same_half_plane(const double *p, const double *q, const double *v1,
                         const double *v2) {
  // Projection on (x,y)-plane
  if (sign_orient2d(p, v1, v2) != sign_orient2d(q, v1, v2))
    return 0;

  // Projection on (y,z)-plane
  if (sign_orient2d(p + 1, v1 + 1, v2 + 1) != sign_orient2d(q + 1, v1 + 1, v2 + 1))
    return 0;

  // Projection on (x,z)-plane
  const double pxz[] = {p[0], p[2]};
  const double qxz[] = {q[0], q[2]};
  const double v1xz[] = {v1[0], v1[2]};
  const double v2xz[] = {v2[0], v2[2]};
  return (sign_orient2d(pxz, v1xz, v2xz) == sign_orient2d(qxz, v1xz, v2xz));
}

//  Input: segments u1-u2 and v1-v2 through their coordinates:
//         u1=(u1x,u1y,u1z), u2=(u2x,u2y,u2z), v1=(v1x,v1y,v1z),
//         v2=(v2x,v2y,v2z).
// Output: 1 -> segments properly intesects, i.e. intersection occours in both
// segments interior,
//         0 -> otherwise.
// Note. Collinear overlapping segments are not considered to be properly
// intersecting.
uint32_t innerSegmentsCross(const double *u1, const double *u2,
                            const double *v1, const double *v2) {

  // Segments have not to share any endpoints
  if (same_point(u1, v1) || same_point(u1, v2) || same_point(u2, v1) ||
      same_point(u2, v2))
    return 0;

  // The 4 endpoints must be coplanar
  if (orient3d(u1, u2, v1, v2) != 0.)
    return 0;

  // Endpoints of one segment cannot stay either on the same side of the other
  // one.
  if (same_half_plane(u1, u2, v1, v2) || same_half_plane(v1, v2, u1, u2))
    return 0;

  // Each segment endpoint cannot be aligned with the other segment.
  if (misAlignment(u1, v1, v2) == 0 || misAlignment(u2, v1, v2) == 0 ||
      misAlignment(v1, u1, u2) == 0 || misAlignment(v2, u1, u2) == 0)
    return 0;

  // If the segment projected on one coordinate plane cross -> segmant cross.
  // Projection on (x,y)-plane
  if (orient2d(u1, u2, v1) != 0. || orient2d(v1, v2, u2) != 0.)
    return 1;

  // Projection on (y,z)-plane
  if (orient2d(u1 + 1, u2 + 1, v1 + 1) != 0. || orient2d(v1 + 1, v2 + 1, u2 + 1) != 0.)
    return 1;

  // Projection on (z,x)-plane
  const double u1xz[] = {u1[0], u1[2]};
  const double u2xz[] = {u2[0], u2[2]};
  const double v1xz[] = {v1[0], v1[2]};
  const double v2xz[] = {v2[0], v2[2]};

  if (orient2d(u1xz, u2xz, v1xz) != 0. || orient2d(v1xz, v2xz, u2xz) != 0.)
    return 1;

  return 0;
}

// ----- Class3: the highest dimensional object is a TRIANGLE -------

//  Input: point p and triangle <v1,v2,v3> through their coordinates:
//         p=(px,py,pz), v1=(v1x,v1y,v1z), v2=(v2x,v2y,v2z), v3=(v3x,v3y,v3z).
// Output: 1 -> point belong to the interior of the triangle,
//         0 -> otherwise.
uint32_t pointInInnerTriangle(const double *p, const double *v1,
                              const double *v2, const double *v3) {
  double o1, o2, oo2, oo4, oo6;

  // Projection on (x,y)-plane -> p VS v1
  o1 = sign_orient2d(p, v2, v3);
  oo2 = o2 = sign_orient2d(v1, v2, v3);
  if (o1 != o2)
    return 0;

  // Projection on (y,z)-plane -> p VS v1
  o1 = sign_orient2d(p + 1, v2 + 1, v3 + 1);
  oo4 = o2 = sign_orient2d(v1 + 1, v2 + 1, v3 + 1);
  if (o1 != o2)
    return 0;

  // Projection on (x,z)-plane -> p VS v1
  const double pxz[] = {p[0], p[2]};
  const double v1xz[] = {v1[0], v1[2]};
  const double v2xz[] = {v2[0], v2[2]};
  const double v3xz[] = {v3[0], v3[2]};
  o1 = sign_orient2d(pxz, v2xz, v3xz);
  oo6 = o2 = sign_orient2d(v1xz, v2xz, v3xz);
  if (o1 != o2)
    return 0;

  // Projection on (x,y)-plane -> p VS v2
  o1 = sign_orient2d(p, v3, v1);
  if (o1 != oo2)
    return 0;

  // Projection on (y,z)-plane -> p VS v2
  o1 = sign_orient2d(p + 1, v3 + 1, v1 + 1);
  if (o1 != oo4)
    return 0;

  // Projection on (x,z)-plane -> p VS v2
  o1 = sign_orient2d(pxz, v3xz, v1xz);
  if (o1 != oo6)
    return 0;

  // Projection on (x,y)-plane -> p VS v3
  o1 = sign_orient2d(p, v1, v2);
  if (o1 != oo2)
    return 0;

  // Projection on (y,z)-plane -> p VS v3
  o1 = sign_orient2d(p + 1, v1 + 1, v2 + 1);
  if (o1 != oo4)
    return 0;

  // Projection on (x,z)-plane -> p VS v3
  o1 = sign_orient2d(pxz, v1xz, v2xz);
  if (o1 != oo6)
    return 0;

  return 1;
}

//  Input: point p and triangle <v1,v2,v3> through their coordinates:
//         p=(px,py,pz), v1=(v1x,v1y,v1z), v2=(v2x,v2y,v2z), v3=(v3x,v3y,v3z).
// Output: 1 -> point belong to the  triangle (interior or boundary),
//         0 -> otherwise.
uint32_t pointInTriangle(const double *p, const double *v1, const double *v2,
                         const double *v3) {

  return (pointInSegment(p, v1, v2) || pointInSegment(p, v2, v3) ||
          pointInSegment(p, v3, v1) || pointInInnerTriangle(p, v1, v2, v3));
}

//  Input: segment u1-u2 and triangle <v1,v2,v3> through their coordinates:
//         u1=(u1x,u1y,u1z), u2=(u2x,u2y,u2z), v1=(v1x,v1y,v1z),
//         v2=(v2x,v2y,v2z), v3=(v3x,v3y,v3z).
// Output: 1 -> segment and triangle properly intesects, i.e. intersection
// occours in both segment and triangle interior,
//        0 -> otherwise.
// Note. Collinear overlapping segments are not considered to be properly
// intersecting.
uint32_t innerSegmentCrossesInnerTriangle(const double *u1, const double *u2,
                                          const double *v1, const double *v2,
                                          const double *v3) {
  // "out of the Box" check.
  double bound;
  bound = std::min(u1[0], u2[0]); // min(u1,u2) alogng x-axis
  if (v1[0] <= bound && v2[0] <= bound && v3[0] <= bound)
    return 0;
  bound = std::max(u1[0], u2[0]); // max(u1,u2) alogng x-axis
  if (v1[0] >= bound && v2[0] >= bound && v3[0] >= bound)
    return 0;
  bound = std::min(u1[1], u2[1]); // min(u1,u2) alogng y-axis
  if (v1[1] <= bound && v2[1] <= bound && v3[1] <= bound)
    return 0;
  bound = std::max(u1[1], u2[1]); // max(u1,u2) alogng y-axis
  if (v1[1] >= bound && v2[1] >= bound && v3[1] >= bound)
    return 0;
  bound = std::min(u1[2], u2[2]); // min(u1,u2) alogng z-axis
  if (v1[2] <= bound && v2[2] <= bound && v3[2] <= bound)
    return 0;
  bound = std::max(u1[2], u2[2]); // max(u1,u2) alogng z-axis
  if (v1[2] >= bound && v2[2] >= bound && v3[2] >= bound)
    return 0;

  const int orient_u1_tri = sign_orient3d(u1, v1, v2, v3);
  const int orient_u2_tri = sign_orient3d(u2, v1, v2, v3);

  // Check if triangle vertices and at least one of the segment endpoints are
  // coplanar: in this case there is no proper intersection.
  if (orient_u1_tri == 0 || orient_u2_tri == 0)
    return 0;

  // Endpoints of one segment cannot stay both in one of the same half-space
  // defined by the triangle.
  if (orient_u1_tri == orient_u2_tri)
    return 0;

  // Since now, endpoints are one abouve and one below the triangle-plane.

  // Intersection between segment and triangle sides are not proper.
  // Check also if segment intersect the triangle-plane outside the triangle.
  const int orient_u_v1v2 = sign_orient3d(u1, u2, v1, v2);
  const int orient_u_v2v3 = sign_orient3d(u1, u2, v2, v3);

  if (orient_u_v1v2 == 0 || orient_u_v2v3 == 0)
    return 0;
  if (orient_u_v1v2 != orient_u_v2v3)
    return 0;

  const int orient_u_v3v1 = sign_orient3d(u1, u2, v3, v1);

  if (orient_u_v3v1 == 0)
    return 0;
  if (orient_u_v3v1 != orient_u_v2v3)
    return 0;

  // Finally, we have a proper intersection.
  return 1;
}

//  Input: segment u1-u2 and triangle <v1,v2,v3> through their coordinates:
//         u1=(u1x,u1y,u1z), u2=(u2x,u2y,u2z),
//         v1=(v1x,v1y,v1z), v2=(v2x,v2y,v2z), v3=(v3x,v3y,v3z).
// Output: 1 -> inner segment and triangle intesects,
//              i.e. intersection occours between segment interior and triangle,
//         0 -> otherwise.
uint32_t innerSegmentCrossesTriangle(const double *u1, const double *u2,
                                     const double *v1, const double *v2,
                                     const double *v3) {
  if (pointInInnerSegment(v1, u1, u2) ||
      pointInInnerSegment(v2, u1, u2) ||
      pointInInnerSegment(v3, u1, u2) ||
      innerSegmentsCross(v2, v3, u1, u2) ||
      innerSegmentsCross(v3, v1, u1, u2) ||
      innerSegmentsCross(v1, v2, u1, u2) ||
      innerSegmentCrossesInnerTriangle(u1, u2, v1, v2, v3))
    return 1;
  
  return 0;
}
