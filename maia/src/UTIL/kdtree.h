// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "functions.h"
#include "pointbox.h"

#ifndef KDTREE_DEFINED
#define KDTREE_DEFINED

/**
 *   Building a K-D-tree
 *     Numercial Recipies in C: The Art of Scientific Computing Thrid Edition
 *     Authors:
 *
 *
 */
template <class T>
inline void SWAP(T& a, T& b) {
  T dum = a;
  a = b;
  b = dum;
}

template <MInt DIM>
struct Boxnode : Box<DIM> {
  MInt mom, dau1, dau2, ptlo, pthi;
  Boxnode() {}
  Boxnode(Point<DIM> mylo, Point<DIM> myhi, MInt mymom, MInt myd1, MInt myd2, MInt myptlo, MInt mypthi)
    : Box<DIM>(mylo, myhi), mom(mymom), dau1(myd1), dau2(myd2), ptlo(myptlo), pthi(mypthi) {}
};
template <MInt DIM>
MInt selecti(const MInt k, MInt* indx, MInt n, MFloat* arr) {
  MInt i, ia, ir, j, l, mid;
  MFloat a;

  l = 0;
  ir = n - 1;
  for(;;) {
    if(ir <= l + 1) {
      if(ir == l + 1 && arr[indx[ir]] < arr[indx[l]]) SWAP(indx[l], indx[ir]);
      return indx[k];
    } else {
      mid = (l + ir) >> 1;
      SWAP(indx[mid], indx[l + 1]);
      if(arr[indx[l]] > arr[indx[ir]]) SWAP(indx[l], indx[ir]);
      if(arr[indx[l + 1]] > arr[indx[ir]]) SWAP(indx[l + 1], indx[ir]);
      if(arr[indx[l]] > arr[indx[l + 1]]) SWAP(indx[l], indx[l + 1]);
      i = l + 1;
      j = ir;
      ia = indx[l + 1];
      a = arr[ia];
      for(;;) {
        do
          i++;
        while(arr[indx[i]] < a);
        do
          j--;
        while(arr[indx[j]] > a);
        if(j < i) break;
        SWAP(indx[i], indx[j]);
      }
      indx[l + 1] = indx[j];
      indx[j] = ia;
      if(j >= k) ir = j - 1;
      if(j <= k) l = i;
    }
  }
}
template <MInt DIM>
struct KDtree {
  static const MFloat BIG;
  MLong nboxes;

  std::vector<Point<DIM>> ptss;
  Boxnode<DIM>* boxes;
  MLong npts;

  std::vector<MInt> ptindx, rptindx;
  //  MInt ccellId;
  KDtree(std::vector<Point<DIM>>& pts);
  ~KDtree() {
    if(boxes != nullptr) delete[] boxes;
    ptss.clear();
  }
  MFloat disti(MInt jpt, MInt kpt);
  MInt locate(Point<DIM> pt);
  MInt locate(MInt jpt);
  MInt nearest(Point<DIM> pt, MFloat& dist);
  void nnearest(MInt jpt, MInt* nn, MFloat* dn, MInt n);
  static void sift_down(MFloat* heap, MInt* ndx, MInt nn);
  MInt locatenear(Point<DIM> pt, MFloat r, MInt* list, MInt nmax, MBool returnCellId = true);
  MInt locatenearest(Point<DIM> pt, MFloat r, MInt* list, MFloat* dn, MInt nmax, MBool& overflow);
};

template <MInt DIM>
const MFloat KDtree<DIM>::BIG(1.0e99);
template <MInt DIM>
KDtree<DIM>::KDtree(std::vector<Point<DIM>>& pts)
  : ptss(pts), boxes(nullptr), npts(pts.size()), ptindx(npts), rptindx(npts) {
  MInt ntmp, m, k, kk, j, nowtask, jbox, np, tmom, tdim, ptlo, pthi;
  MInt* hp;
  MFloat* cp;
  MInt taskmom[50], taskdim[50];
  if(npts == 0) return;
  for(k = 0; k < npts; k++)
    ptindx[k] = k;
  m = 1;
  for(ntmp = npts; ntmp; ntmp >>= 1) {
    m <<= 1;
  }
  nboxes = 2 * npts - (m >> 1);
  if(m < nboxes) nboxes = m;
  nboxes--;
  if(nboxes > 0) {
    boxes = new Boxnode<DIM>[nboxes];
    MFloat* coords = new MFloat[DIM * npts];
    for(j = 0, kk = 0; j < DIM; j++, kk += npts) {
      for(k = 0; k < npts; k++)
        coords[kk + k] = pts[k].x[j];
    }
    Point<DIM> lo(-BIG, -BIG, -BIG), hi(BIG, BIG, BIG);
    boxes[0] = Boxnode<DIM>(lo, hi, 0, 0, 0, 0, npts - 1);
    jbox = 0;
    taskmom[1] = 0;
    taskdim[1] = 0;
    nowtask = 1;
    while(nowtask) {
      tmom = taskmom[nowtask];
      tdim = taskdim[nowtask--];
      ptlo = boxes[tmom].ptlo;
      pthi = boxes[tmom].pthi;
      hp = &ptindx[ptlo];
      cp = &coords[tdim * npts];
      np = pthi - ptlo + 1;
      kk = (np - 1) / 2;
      (void)selecti<DIM>(kk, hp, np, cp);
      hi = boxes[tmom].hi;
      lo = boxes[tmom].lo;
      hi.x[tdim] = lo.x[tdim] = coords[tdim * npts + hp[kk]];
      if(jbox < nboxes - 1) boxes[++jbox] = Boxnode<DIM>(boxes[tmom].lo, hi, tmom, 0, 0, ptlo, ptlo + kk);
      if(jbox < nboxes - 1) boxes[++jbox] = Boxnode<DIM>(lo, boxes[tmom].hi, tmom, 0, 0, ptlo + kk + 1, pthi);
      if(tmom < nboxes) boxes[tmom].dau1 = jbox - 1;
      if(tmom < nboxes) boxes[tmom].dau2 = jbox;
      if(kk > 1) {
        taskmom[++nowtask] = jbox - 1;
        taskdim[nowtask] = (tdim + 1) % DIM;
      }
      if(np - kk > 3) {
        taskmom[++nowtask] = jbox;
        taskdim[nowtask] = (tdim + 1) % DIM;
      }
    }
    for(j = 0; j < npts; j++)
      rptindx[ptindx[j]] = j;
    delete[] coords;
  }
}
template <MInt DIM>
MFloat KDtree<DIM>::disti(MInt jpt, MInt kpt) {
  if(jpt == kpt)
    return BIG;
  else
    return dist(ptss[jpt], ptss[kpt]);
}
template <MInt DIM>
MInt KDtree<DIM>::locate(Point<DIM> pt) {
  MInt nb, d1, jdim;
  nb = jdim = 0;
  while(boxes[nb].dau1) {
    d1 = boxes[nb].dau1;
    if(pt.x[jdim] <= boxes[d1].hi.x[jdim])
      nb = d1;
    else
      nb = boxes[nb].dau2;
    //		jdim = ++jdim % DIM;
    ++jdim;
    jdim %= DIM;
  }
  return nb;
}
template <MInt DIM>
MInt KDtree<DIM>::locate(MInt jpt) {
  MInt nb, d1, jh;
  jh = rptindx[jpt];
  nb = 0;
  while(boxes[nb].dau1) {
    d1 = boxes[nb].dau1;
    if(jh <= boxes[d1].pthi)
      nb = d1;
    else
      nb = boxes[nb].dau2;
  }
  return nb;
}
template <MInt DIM>
MInt KDtree<DIM>::nearest(Point<DIM> pt, MFloat& dist1) {
  MInt i, k, ntask;
  MInt nrst = 0;
  MInt task[50];
  MFloat dnrst = BIG, d;
  k = locate(pt);
  for(i = boxes[k].ptlo; i <= boxes[k].pthi; i++) {
    d = dist(ptss[ptindx[i]], pt);
    if(d < dnrst) {
      nrst = ptindx[i];
      dnrst = d;
    }
  }
  task[1] = 0;
  ntask = 1;
  while(ntask) {
    k = task[ntask--];
    if(dist(boxes[k], pt) < dnrst) {
      if(boxes[k].dau1) {
        task[++ntask] = boxes[k].dau1;
        task[++ntask] = boxes[k].dau2;
      } else {
        for(i = boxes[k].ptlo; i <= boxes[k].pthi; i++) {
          d = dist(ptss[ptindx[i]], pt);
          if(d < dnrst) {
            nrst = ptindx[i];
            dnrst = d;
          }
        }
      }
    }
  }
  dist1 = dist(ptss[nrst], pt); // gives the distance!
  return ptss[nrst].cellId;     // nrst;
}
template <MInt DIM>
void KDtree<DIM>::nnearest(MInt jpt, MInt* nn, MFloat* dn, MInt n) {
  MInt i, k, ntask, kp;
  MInt task[50];
  MFloat d;
  ASSERT(!(n > npts - 1), "too many neighbors requested");
  for(i = 0; i < n; i++)
    nn[i] = -1;
  for(i = 0; i < n; i++)
    dn[i] = BIG;
  if(npts == 0) {
    return;
  }
  if(npts == 1) {
    nn[0] = ptss[ptindx[0]].cellId;
    dn[0] = disti(ptindx[0], jpt);
    return;
  } else if(npts == 2) {
    MInt smaller = (disti(ptindx[0], jpt) < disti(ptindx[1], jpt)) ? 0 : 1;
    MInt other = (smaller + 1) % 2;
    nn[0] = ptss[ptindx[smaller]].cellId;
    dn[0] = disti(ptindx[smaller], jpt);
    if(n > 1) {
      nn[1] = ptss[ptindx[other]].cellId;
      dn[1] = disti(ptindx[other], jpt);
    }
    return;
  }
  kp = boxes[locate(jpt)].mom;
  while(boxes[kp].pthi - boxes[kp].ptlo < n)
    kp = boxes[kp].mom;
  for(i = boxes[kp].ptlo; i <= boxes[kp].pthi; i++) {
    if(jpt == ptindx[i]) continue;
    d = disti(ptindx[i], jpt);
    if(d < dn[0]) {
      dn[0] = d;
      nn[0] = ptindx[i];
      if(n > 1) sift_down(dn, nn, n);
    }
  }
  task[1] = 0;
  ntask = 1;
  while(ntask) {
    k = task[ntask--];
    if(k == kp) continue;
    if(dist(boxes[k], ptss[jpt]) < dn[0]) {
      if(boxes[k].dau1) {
        task[++ntask] = boxes[k].dau1;
        task[++ntask] = boxes[k].dau2;
      } else {
        for(i = boxes[k].ptlo; i <= boxes[k].pthi; i++) {
          d = disti(ptindx[i], jpt);
          if(d < dn[0]) {
            dn[0] = d;
            nn[0] = ptindx[i];
            if(n > 1) sift_down(dn, nn, n);
          }
        }
      }
    }
  }
  for(i = 0; i < n; i++)
    nn[i] = (nn[i] > -1) ? ptss[nn[i]].cellId : -1;
  return;
}
template <MInt DIM>
void KDtree<DIM>::sift_down(MFloat* heap, MInt* ndx, MInt nn) {
  MInt n = nn - 1;
  MInt j, jold, ia;
  MFloat a;
  a = heap[0];
  ia = ndx[0];
  jold = 0;
  j = 1;
  while(j <= n) {
    if(j < n && heap[j] < heap[j + 1]) j++;
    if(a >= heap[j]) break;
    heap[jold] = heap[j];
    ndx[jold] = ndx[j];
    jold = j;
    j = 2 * j + 1;
  }
  heap[jold] = a;
  ndx[jold] = ia;
}
template <MInt DIM>
MInt KDtree<DIM>::locatenear(Point<DIM> pt, MFloat r, MInt* list, MInt nmax, MBool returnCellId) {
  MInt k, i, nb, nbold, nret, ntask, jdim, d1, d2;
  MInt task[50];
  nb = jdim = nret = 0;
  ASSERT(!(r < 0.0), "radius must be nonnegative");
  ASSERT(nmax > 0, "");
  if(npts == 0) {
    return 0;
  } else if(npts == 1) {
    if(dist(ptss[ptindx[0]], pt) < r && nmax > 0) list[nret++] = ptss[ptindx[0]].cellId;
  } else if(npts == 2) {
    if(dist(ptss[ptindx[0]], pt) < r && nmax > 0) list[nret++] = ptss[ptindx[0]].cellId;
    if(dist(ptss[ptindx[1]], pt) < r && nmax > 1) list[nret++] = ptss[ptindx[1]].cellId;
  } else {
    while(boxes[nb].dau1) {
      nbold = nb;
      d1 = boxes[nb].dau1;
      d2 = boxes[nb].dau2;
      if(pt.x[jdim] + r <= boxes[d1].hi.x[jdim])
        nb = d1;
      else if(pt.x[jdim] - r >= boxes[d2].lo.x[jdim])
        nb = d2;
      // jdim = ++jdim % DIM;//undefined behavior
      jdim = (jdim + 1) % DIM;
      if(nb == nbold) break;
    }
    task[1] = nb;
    ntask = 1;
    while(ntask) {
      k = task[ntask--];
      if(dist(boxes[k], pt) > r) continue;
      if(boxes[k].dau1) {
        task[++ntask] = boxes[k].dau1;
        task[++ntask] = boxes[k].dau2;
      } else {
        for(i = boxes[k].ptlo; i <= boxes[k].pthi; i++) {
          if(dist(ptss[ptindx[i]], pt) <= r && nret < nmax) {
            if(returnCellId) {
              list[nret++] = ptss[ptindx[i]].cellId;
            } else {
              list[nret++] = ptindx[i];
            }
          }

          if(nret == nmax) return nmax;
        }
      }
    }
  }
  return nret;
}
template <MInt DIM>
MInt KDtree<DIM>::locatenearest(Point<DIM> pt, MFloat r, MInt* list, MFloat* dn, MInt nmax, MBool& overflow) {
  MInt k, i, nb, nbold, nret, ntask, jdim, d1, d2;
  MInt task[50];
  nb = jdim = nret = 0;
  overflow = false;
  ASSERT(!(r < 0.0), "radius must be nonnegative");
  ASSERT(nmax > 0, "");
  if(npts == 0) {
    return 0;
  } else if(npts == 1) {
    if(dist(ptss[ptindx[0]], pt) < r && nmax > 0) list[nret++] = ptss[ptindx[0]].cellId;
  } else if(npts == 2) {
    MInt smaller = (dist(ptss[ptindx[0]], pt) < dist(ptss[ptindx[1]], pt)) ? 0 : 1;
    MInt other = (smaller + 1) % 2;
    if(dist(ptss[ptindx[smaller]], pt) < r && nmax > 0) list[nret++] = ptss[ptindx[smaller]].cellId;
    if(dist(ptss[ptindx[other]], pt) < r && nmax > 1) list[nret++] = ptss[ptindx[other]].cellId;
  } else {
    while(boxes[nb].dau1) {
      nbold = nb;
      d1 = boxes[nb].dau1;
      d2 = boxes[nb].dau2;
      if(pt.x[jdim] + r <= boxes[d1].hi.x[jdim])
        nb = d1;
      else if(pt.x[jdim] - r >= boxes[d2].lo.x[jdim])
        nb = d2;
      // jdim = ++jdim % DIM;//undefined behavior
      jdim = (jdim + 1) % DIM;
      if(nb == nbold) break;
    }
    task[1] = nb;
    ntask = 1;
    while(ntask) {
      k = task[ntask--];
      if(dist(boxes[k], pt) > r) continue;
      if(boxes[k].dau1) {
        task[++ntask] = boxes[k].dau1;
        task[++ntask] = boxes[k].dau2;
      } else {
        for(i = boxes[k].ptlo; i <= boxes[k].pthi; i++) {
          MFloat dst = dist(ptss[ptindx[i]], pt);
          if(dst < r) {
            if(nret == nmax) overflow = true;
            MInt s = 0;
            for(s = 0; s < nret; s++) {
              if(dst < dn[s]) {
                for(MInt t = std::min(nret, nmax - 1); t > s; t--) {
                  list[t] = list[t - 1];
                  dn[t] = dn[t - 1];
                }
                break;
              }
            }
            if(s < nmax) {
              list[s] = ptss[ptindx[i]].cellId;
              dn[s] = dst;
            }
            if(nret < nmax) nret++;
          }
        }
      }
    }
  }
  return nret;
}

#endif
