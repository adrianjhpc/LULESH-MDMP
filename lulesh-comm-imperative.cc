#include "lulesh.h"

// If no MPI, then this whole file is stubbed out
#if USE_MPI

#include <string.h>

#include "mdmp_interface.h"

// Note: Explicit mdmp_wait declaration removed. 
// The LLVM compiler pass handles wait injection automatically now!

/* Comm Routines */

#define ALLOW_UNPACKED_PLANE false
#define ALLOW_UNPACKED_ROW   false
#define ALLOW_UNPACKED_COL   false

/******************************************/

void CommRecv(Domain& domain, Int_t msgType, Index_t xferFields,
              Index_t dx, Index_t dy, Index_t dz, bool doRecv, bool planeOnly) {

   if (domain.numRanks() == 1)
      return ;

   int myRank = MDMP_GET_RANK() ;
   Index_t maxPlaneComm = xferFields * domain.maxPlaneSize() ;
   Index_t maxEdgeComm  = xferFields * domain.maxEdgeSize() ;
   Index_t pmsg = 0 ; /* plane comm msg */
   Index_t emsg = 0 ; /* edge comm msg */
   Index_t cmsg = 0 ; /* corner comm msg */
   bool rowMin, rowMax, colMin, colMax, planeMin, planeMax ;

   /* assume communication to 6 neighbors by default */
   rowMin = rowMax = colMin = colMax = planeMin = planeMax = true ;

   if (domain.rowLoc() == 0) rowMin = false ;
   if (domain.rowLoc() == (domain.tp()-1)) rowMax = false ;
   if (domain.colLoc() == 0) colMin = false ;
   if (domain.colLoc() == (domain.tp()-1)) colMax = false ;
   if (domain.planeLoc() == 0) planeMin = false ;
   if (domain.planeLoc() == (domain.tp()-1)) planeMax = false ;

   if (planeMin && doRecv) {
      int fromRank = myRank - domain.tp()*domain.tp() ;
      int recvCount = dx * dy * xferFields ;
      MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm], recvCount, myRank, fromRank, msgType);
      ++pmsg ;
   }
   if (planeMax) {
      int fromRank = myRank + domain.tp()*domain.tp() ;
      int recvCount = dx * dy * xferFields ;
      MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm], recvCount, myRank, fromRank, msgType);
      ++pmsg ;
   }
   if (rowMin && doRecv) {
      int fromRank = myRank - domain.tp() ;
      int recvCount = dx * dz * xferFields ;
      MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm], recvCount, myRank, fromRank, msgType);
      ++pmsg ;
   }
   if (rowMax) {
      int fromRank = myRank + domain.tp() ;
      int recvCount = dx * dz * xferFields ;
      MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm], recvCount, myRank, fromRank, msgType);
      ++pmsg ;
   }
   if (colMin && doRecv) {
      int fromRank = myRank - 1 ;
      int recvCount = dy * dz * xferFields ;
      MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm], recvCount, myRank, fromRank, msgType);
      ++pmsg ;
   }
   if (colMax) {
      int fromRank = myRank + 1 ;
      int recvCount = dy * dz * xferFields ;
      MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm], recvCount, myRank, fromRank, msgType);
      ++pmsg ;
   }

   if (!planeOnly) {
      /* receive data from domains connected only by an edge */
      if (rowMin && colMin && doRecv) {
         int fromRank = myRank - domain.tp() - 1 ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm], dz * xferFields, myRank, fromRank, msgType);
         ++emsg ;
      }
      if (rowMin && planeMin && doRecv) {
         int fromRank = myRank - domain.tp()*domain.tp() - domain.tp() ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm], dx * xferFields, myRank, fromRank, msgType);
         ++emsg ;
      }
      if (colMin && planeMin && doRecv) {
         int fromRank = myRank - domain.tp()*domain.tp() - 1 ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm], dy * xferFields, myRank, fromRank, msgType);
         ++emsg ;
      }
      if (rowMax && colMax) {
         int fromRank = myRank + domain.tp() + 1 ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm], dz * xferFields, myRank, fromRank, msgType);
         ++emsg ;
      }
      if (rowMax && planeMax) {
         int fromRank = myRank + domain.tp()*domain.tp() + domain.tp() ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm], dx * xferFields, myRank, fromRank, msgType);
         ++emsg ;
      }
      if (colMax && planeMax) {
         int fromRank = myRank + domain.tp()*domain.tp() + 1 ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm], dy * xferFields, myRank, fromRank, msgType);
         ++emsg ;
      }
      if (rowMax && colMin) {
         int fromRank = myRank + domain.tp() - 1 ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm], dz * xferFields, myRank, fromRank, msgType);
         ++emsg ;
      }
      if (rowMin && planeMax) {
         int fromRank = myRank + domain.tp()*domain.tp() - domain.tp() ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm], dx * xferFields, myRank, fromRank, msgType);
         ++emsg ;
      }
      if (colMin && planeMax) {
         int fromRank = myRank + domain.tp()*domain.tp() - 1 ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm], dy * xferFields, myRank, fromRank, msgType);
         ++emsg ;
      }
      if (rowMin && colMax && doRecv) {
         int fromRank = myRank - domain.tp() + 1 ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm], dz * xferFields, myRank, fromRank, msgType);
         ++emsg ;
      }
      if (rowMax && planeMin && doRecv) {
         int fromRank = myRank - domain.tp()*domain.tp() + domain.tp() ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm], dx * xferFields, myRank, fromRank, msgType);
         ++emsg ;
      }
      if (colMax && planeMin && doRecv) {
         int fromRank = myRank - domain.tp()*domain.tp() + 1 ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm], dy * xferFields, myRank, fromRank, msgType);
         ++emsg ;
      }

      /* receive data from domains connected only by a corner */
      if (rowMin && colMin && planeMin && doRecv) {
         int fromRank = myRank - domain.tp()*domain.tp() - domain.tp() - 1 ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL], xferFields, myRank, fromRank, msgType);
         ++cmsg ;
      }
      if (rowMin && colMin && planeMax) {
         int fromRank = myRank + domain.tp()*domain.tp() - domain.tp() - 1 ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL], xferFields, myRank, fromRank, msgType);
         ++cmsg ;
      }
      if (rowMin && colMax && planeMin && doRecv) {
         int fromRank = myRank - domain.tp()*domain.tp() - domain.tp() + 1 ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL], xferFields, myRank, fromRank, msgType);
         ++cmsg ;
      }
      if (rowMin && colMax && planeMax) {
         int fromRank = myRank + domain.tp()*domain.tp() - domain.tp() + 1 ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL], xferFields, myRank, fromRank, msgType);
         ++cmsg ;
      }
      if (rowMax && colMin && planeMin && doRecv) {
         int fromRank = myRank - domain.tp()*domain.tp() + domain.tp() - 1 ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL], xferFields, myRank, fromRank, msgType);
         ++cmsg ;
      }
      if (rowMax && colMin && planeMax) {
         int fromRank = myRank + domain.tp()*domain.tp() + domain.tp() - 1 ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL], xferFields, myRank, fromRank, msgType);
         ++cmsg ;
      }
      if (rowMax && colMax && planeMin && doRecv) {
         int fromRank = myRank - domain.tp()*domain.tp() + domain.tp() + 1 ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL], xferFields, myRank, fromRank, msgType);
         ++cmsg ;
      }
      if (rowMax && colMax && planeMax) {
         int fromRank = myRank + domain.tp()*domain.tp() + domain.tp() + 1 ;
         MDMP_RECV(&domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL], xferFields, myRank, fromRank, msgType);
         ++cmsg ;
      }
   }
}

/******************************************/

void CommSend(Domain& domain, Int_t msgType,
              Index_t xferFields, Domain_member *fieldData,
              Index_t dx, Index_t dy, Index_t dz, bool doSend, bool planeOnly)
{
   if (domain.numRanks() == 1)
      return ;

   int myRank = MDMP_GET_RANK() ;
   Index_t maxPlaneComm = xferFields * domain.maxPlaneSize() ;
   Index_t maxEdgeComm  = xferFields * domain.maxEdgeSize() ;
   Index_t pmsg = 0 ; /* plane comm msg */
   Index_t emsg = 0 ; /* edge comm msg */
   Index_t cmsg = 0 ; /* corner comm msg */
   Real_t *destAddr ;
   bool rowMin, rowMax, colMin, colMax, planeMin, planeMax ;
   
   rowMin = rowMax = colMin = colMax = planeMin = planeMax = true ;
   if (domain.rowLoc() == 0) rowMin = false ;
   if (domain.rowLoc() == (domain.tp()-1)) rowMax = false ;
   if (domain.colLoc() == 0) colMin = false ;
   if (domain.colLoc() == (domain.tp()-1)) colMax = false ;
   if (domain.planeLoc() == 0) planeMin = false ;
   if (domain.planeLoc() == (domain.tp()-1)) planeMax = false ;

   if (planeMin | planeMax) {
      int sendCount = dx * dy ;
      if (planeMin) {
         destAddr = &domain.commDataSend[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi] ;
            for (Index_t i=0; i<sendCount; ++i) {
               destAddr[i] = (domain.*src)(i) ;
            }
            destAddr += sendCount ;
         }
         destAddr -= xferFields*sendCount ;
         MDMP_SEND(destAddr, xferFields*sendCount, myRank, myRank - domain.tp()*domain.tp(), msgType);
         ++pmsg ;
      }
      if (planeMax && doSend) {
         destAddr = &domain.commDataSend[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi] ;
            for (Index_t i=0; i<sendCount; ++i) {
               destAddr[i] = (domain.*src)(dx*dy*(dz - 1) + i) ;
            }
            destAddr += sendCount ;
         }
         destAddr -= xferFields*sendCount ;
         MDMP_SEND(destAddr, xferFields*sendCount, myRank, myRank + domain.tp()*domain.tp(), msgType);
         ++pmsg ;
      }
   }
   if (rowMin | rowMax) {
      int sendCount = dx * dz ;
      if (rowMin) {
         destAddr = &domain.commDataSend[pmsg * maxPlaneComm] ;
         for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi] ;
            for (Index_t i=0; i<dz; ++i) {
               for (Index_t j=0; j<dx; ++j) {
                  destAddr[i*dx+j] = (domain.*src)(i*dx*dy + j) ;
               }
            }
            destAddr += sendCount ;
         }
         destAddr -= xferFields*sendCount ;
         MDMP_SEND(destAddr, xferFields*sendCount, myRank, myRank - domain.tp(), msgType);
         ++pmsg ;
      }
      if (rowMax && doSend) {
         destAddr = &domain.commDataSend[pmsg * maxPlaneComm] ;
         for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi] ;
            for (Index_t i=0; i<dz; ++i) {
               for (Index_t j=0; j<dx; ++j) {
                  destAddr[i*dx+j] = (domain.*src)(dx*(dy - 1) + i*dx*dy + j) ;
               }
            }
            destAddr += sendCount ;
         }
         destAddr -= xferFields*sendCount ;
         MDMP_SEND(destAddr, xferFields*sendCount, myRank, myRank + domain.tp(), msgType);
         ++pmsg ;
      }
   }
   if (colMin | colMax) {
      int sendCount = dy * dz ;
      if (colMin) {
         destAddr = &domain.commDataSend[pmsg * maxPlaneComm] ;
         for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi] ;
            for (Index_t i=0; i<dz; ++i) {
               for (Index_t j=0; j<dy; ++j) {
                  destAddr[i*dy + j] = (domain.*src)(i*dx*dy + j*dx) ;
               }
            }
            destAddr += sendCount ;
         }
         destAddr -= xferFields*sendCount ;
         MDMP_SEND(destAddr, xferFields*sendCount, myRank, myRank - 1, msgType);
         ++pmsg ;
      }
      if (colMax && doSend) {
         destAddr = &domain.commDataSend[pmsg * maxPlaneComm] ;
         for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi] ;
            for (Index_t i=0; i<dz; ++i) {
               for (Index_t j=0; j<dy; ++j) {
                  destAddr[i*dy + j] = (domain.*src)(dx - 1 + i*dx*dy + j*dx) ;
               }
            }
            destAddr += sendCount ;
         }
         destAddr -= xferFields*sendCount ;
         MDMP_SEND(destAddr, xferFields*sendCount, myRank, myRank + 1, msgType);
         ++pmsg ;
      }
   }

   if (!planeOnly) {
      if (rowMin && colMin) {
         int toRank = myRank - domain.tp() - 1 ;
         destAddr = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
         for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi] ;
            for (Index_t i=0; i<dz; ++i) {
               destAddr[i] = (domain.*src)(i*dx*dy) ;
            }
            destAddr += dz ;
         }
         destAddr -= xferFields*dz ;
         MDMP_SEND(destAddr, xferFields*dz, myRank, toRank, msgType);
         ++emsg ;
      }
      if (rowMin && planeMin) {
         int toRank = myRank - domain.tp()*domain.tp() - domain.tp() ;
         destAddr = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
         for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi] ;
            for (Index_t i=0; i<dx; ++i) {
               destAddr[i] = (domain.*src)(i) ;
            }
            destAddr += dx ;
         }
         destAddr -= xferFields*dx ;
         MDMP_SEND(destAddr, xferFields*dx, myRank, toRank, msgType);
         ++emsg ;
      }
      if (colMin && planeMin) {
         int toRank = myRank - domain.tp()*domain.tp() - 1 ;
         destAddr = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
         for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi] ;
            for (Index_t i=0; i<dy; ++i) {
               destAddr[i] = (domain.*src)(i*dx) ;
            }
            destAddr += dy ;
         }
         destAddr -= xferFields*dy ;
         MDMP_SEND(destAddr, xferFields*dy, myRank, toRank, msgType);
         ++emsg ;
      }
      if (rowMax && colMax && doSend) {
         int toRank = myRank + domain.tp() + 1 ;
         destAddr = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
         for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi] ;
            for (Index_t i=0; i<dz; ++i) {
               destAddr[i] = (domain.*src)(dx*dy - 1 + i*dx*dy) ;
            }
            destAddr += dz ;
         }
         destAddr -= xferFields*dz ;
         MDMP_SEND(destAddr, xferFields*dz, myRank, toRank, msgType);
         ++emsg ;
      }
      if (rowMax && planeMax && doSend) {
         int toRank = myRank + domain.tp()*domain.tp() + domain.tp() ;
         destAddr = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
         for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi] ;
            for (Index_t i=0; i<dx; ++i) {
              destAddr[i] = (domain.*src)(dx*(dy-1) + dx*dy*(dz-1) + i) ;
            }
            destAddr += dx ;
         }
         destAddr -= xferFields*dx ;
         MDMP_SEND(destAddr, xferFields*dx, myRank, toRank, msgType);
         ++emsg ;
      }
      if (colMax && planeMax && doSend) {
         int toRank = myRank + domain.tp()*domain.tp() + 1 ;
         destAddr = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
         for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi] ;
            for (Index_t i=0; i<dy; ++i) {
               destAddr[i] = (domain.*src)(dx*dy*(dz-1) + dx - 1 + i*dx) ;
            }
            destAddr += dy ;
         }
         destAddr -= xferFields*dy ;
         MDMP_SEND(destAddr, xferFields*dy, myRank, toRank, msgType);
         ++emsg ;
      }
      if (rowMax && colMin && doSend) {
         int toRank = myRank + domain.tp() - 1 ;
         destAddr = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
         for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi] ;
            for (Index_t i=0; i<dz; ++i) {
               destAddr[i] = (domain.*src)(dx*(dy-1) + i*dx*dy) ;
            }
            destAddr += dz ;
         }
         destAddr -= xferFields*dz ;
         MDMP_SEND(destAddr, xferFields*dz, myRank, toRank, msgType);
         ++emsg ;
      }
      if (rowMin && planeMax && doSend) {
         int toRank = myRank + domain.tp()*domain.tp() - domain.tp() ;
         destAddr = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
         for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi] ;
            for (Index_t i=0; i<dx; ++i) {
               destAddr[i] = (domain.*src)(dx*dy*(dz-1) + i) ;
            }
            destAddr += dx ;
         }
         destAddr -= xferFields*dx ;
         MDMP_SEND(destAddr, xferFields*dx, myRank, toRank, msgType);
         ++emsg ;
      }
      if (colMin && planeMax && doSend) {
         int toRank = myRank + domain.tp()*domain.tp() - 1 ;
         destAddr = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
         for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi] ;
            for (Index_t i=0; i<dy; ++i) {
               destAddr[i] = (domain.*src)(dx*dy*(dz-1) + i*dx) ;
            }
            destAddr += dy ;
         }
         destAddr -= xferFields*dy ;
         MDMP_SEND(destAddr, xferFields*dy, myRank, toRank, msgType);
         ++emsg ;
      }
      if (rowMin && colMax) {
         int toRank = myRank - domain.tp() + 1 ;
         destAddr = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
         for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi] ;
            for (Index_t i=0; i<dz; ++i) {
               destAddr[i] = (domain.*src)(dx - 1 + i*dx*dy) ;
            }
            destAddr += dz ;
         }
         destAddr -= xferFields*dz ;
         MDMP_SEND(destAddr, xferFields*dz, myRank, toRank, msgType);
         ++emsg ;
      }
      if (rowMax && planeMin) {
         int toRank = myRank - domain.tp()*domain.tp() + domain.tp() ;
         destAddr = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
         for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi] ;
            for (Index_t i=0; i<dx; ++i) {
               destAddr[i] = (domain.*src)(dx*(dy - 1) + i) ;
            }
            destAddr += dx ;
         }
         destAddr -= xferFields*dx ;
         MDMP_SEND(destAddr, xferFields*dx, myRank, toRank, msgType);
         ++emsg ;
      }
      if (colMax && planeMin) {
         int toRank = myRank - domain.tp()*domain.tp() + 1 ;
         destAddr = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
         for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi] ;
            for (Index_t i=0; i<dy; ++i) {
               destAddr[i] = (domain.*src)(dx - 1 + i*dx) ;
            }
            destAddr += dy ;
         }
         destAddr -= xferFields*dy ;
         MDMP_SEND(destAddr, xferFields*dy, myRank, toRank, msgType);
         ++emsg ;
      }

      if (rowMin && colMin && planeMin) {
         int toRank = myRank - domain.tp()*domain.tp() - domain.tp() - 1 ;
         Real_t *comBuf = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
         for (Index_t fi=0; fi<xferFields; ++fi) comBuf[fi] = (domain.*fieldData[fi])(0) ;
         MDMP_SEND(comBuf, xferFields, myRank, toRank, msgType);
         ++cmsg ;
      }
      if (rowMin && colMin && planeMax && doSend) {
         int toRank = myRank + domain.tp()*domain.tp() - domain.tp() - 1 ;
         Real_t *comBuf = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
         Index_t idx = dx*dy*(dz - 1) ;
         for (Index_t fi=0; fi<xferFields; ++fi) comBuf[fi] = (domain.*fieldData[fi])(idx) ;
         MDMP_SEND(comBuf, xferFields, myRank, toRank, msgType);
         ++cmsg ;
      }
      if (rowMin && colMax && planeMin) {
         int toRank = myRank - domain.tp()*domain.tp() - domain.tp() + 1 ;
         Real_t *comBuf = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
         Index_t idx = dx - 1 ;
         for (Index_t fi=0; fi<xferFields; ++fi) comBuf[fi] = (domain.*fieldData[fi])(idx) ;
         MDMP_SEND(comBuf, xferFields, myRank, toRank, msgType);
         ++cmsg ;
      }
      if (rowMin && colMax && planeMax && doSend) {
         int toRank = myRank + domain.tp()*domain.tp() - domain.tp() + 1 ;
         Real_t *comBuf = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
         Index_t idx = dx*dy*(dz - 1) + (dx - 1) ;
         for (Index_t fi=0; fi<xferFields; ++fi) comBuf[fi] = (domain.*fieldData[fi])(idx) ;
         MDMP_SEND(comBuf, xferFields, myRank, toRank, msgType);
         ++cmsg ;
      }
      if (rowMax && colMin && planeMin) {
         int toRank = myRank - domain.tp()*domain.tp() + domain.tp() - 1 ;
         Real_t *comBuf = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
         Index_t idx = dx*(dy - 1) ;
         for (Index_t fi=0; fi<xferFields; ++fi) comBuf[fi] = (domain.*fieldData[fi])(idx) ;
         MDMP_SEND(comBuf, xferFields, myRank, toRank, msgType);
         ++cmsg ;
      }
      if (rowMax && colMin && planeMax && doSend) {
         int toRank = myRank + domain.tp()*domain.tp() + domain.tp() - 1 ;
         Real_t *comBuf = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
         Index_t idx = dx*dy*(dz - 1) + dx*(dy - 1) ;
         for (Index_t fi=0; fi<xferFields; ++fi) comBuf[fi] = (domain.*fieldData[fi])(idx) ;
         MDMP_SEND(comBuf, xferFields, myRank, toRank, msgType);
         ++cmsg ;
      }
      if (rowMax && colMax && planeMin) {
         int toRank = myRank - domain.tp()*domain.tp() + domain.tp() + 1 ;
         Real_t *comBuf = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
         Index_t idx = dx*dy - 1 ;
         for (Index_t fi=0; fi<xferFields; ++fi) comBuf[fi] = (domain.*fieldData[fi])(idx) ;
         MDMP_SEND(comBuf, xferFields, myRank, toRank, msgType);
         ++cmsg ;
      }
      if (rowMax && colMax && planeMax && doSend) {
         int toRank = myRank + domain.tp()*domain.tp() + domain.tp() + 1 ;
         Real_t *comBuf = &domain.commDataSend[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
         Index_t idx = dx*dy*dz - 1 ;
         for (Index_t fi=0; fi<xferFields; ++fi) comBuf[fi] = (domain.*fieldData[fi])(idx) ;
         MDMP_SEND(comBuf, xferFields, myRank, toRank, msgType);
         ++cmsg ;
      }
   }

   // MDMP_COMMIT(); removed - Imperative sends fire immediately
}

/******************************************/

void CommSBN(Domain& domain, Int_t xferFields, Domain_member *fieldData) {

   if (domain.numRanks() == 1)
      return ;

   // mdmp_wait(-1); removed - LLVM pass intercepts memory read

   Index_t maxPlaneComm = xferFields * domain.maxPlaneSize() ;
   Index_t maxEdgeComm  = xferFields * domain.maxEdgeSize() ;
   Index_t pmsg = 0 ; /* plane comm msg */
   Index_t emsg = 0 ; /* edge comm msg */
   Index_t cmsg = 0 ; /* corner comm msg */
   Index_t dx = domain.sizeX() + 1 ;
   Index_t dy = domain.sizeY() + 1 ;
   Index_t dz = domain.sizeZ() + 1 ;
   Real_t *srcAddr ;
   Index_t rowMin, rowMax, colMin, colMax, planeMin, planeMax ;
   
   rowMin = rowMax = colMin = colMax = planeMin = planeMax = 1 ;
   if (domain.rowLoc() == 0) rowMin = 0 ;
   if (domain.rowLoc() == (domain.tp()-1)) rowMax = 0 ;
   if (domain.colLoc() == 0) colMin = 0 ;
   if (domain.colLoc() == (domain.tp()-1)) colMax = 0 ;
   if (domain.planeLoc() == 0) planeMin = 0 ;
   if (domain.planeLoc() == (domain.tp()-1)) planeMax = 0 ;

   if (planeMin | planeMax) {
      Index_t opCount = dx * dy ;

      if (planeMin) {
         srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi] ;
            for (Index_t i=0; i<opCount; ++i) {
               (domain.*dest)(i) += srcAddr[i] ;
            }
            srcAddr += opCount ;
         }
         ++pmsg ;
      }
      if (planeMax) {
         srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi] ;
            for (Index_t i=0; i<opCount; ++i) {
               (domain.*dest)(dx*dy*(dz - 1) + i) += srcAddr[i] ;
            }
            srcAddr += opCount ;
         }
         ++pmsg ;
      }
   }

   if (rowMin | rowMax) {
      Index_t opCount = dx * dz ;

      if (rowMin) {
         srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi] ;
            for (Index_t i=0; i<dz; ++i) {
               for (Index_t j=0; j<dx; ++j) {
                  (domain.*dest)(i*dx*dy + j) += srcAddr[i*dx + j] ;
               }
            }
            srcAddr += opCount ;
         }
         ++pmsg ;
      }
      if (rowMax) {
         srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi] ;
            for (Index_t i=0; i<dz; ++i) {
               for (Index_t j=0; j<dx; ++j) {
                  (domain.*dest)(dx*(dy - 1) + i*dx*dy + j) += srcAddr[i*dx + j] ;
               }
            }
            srcAddr += opCount ;
         }
         ++pmsg ;
      }
   }
   if (colMin | colMax) {
      Index_t opCount = dy * dz ;

      if (colMin) {
         srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi] ;
            for (Index_t i=0; i<dz; ++i) {
               for (Index_t j=0; j<dy; ++j) {
                  (domain.*dest)(i*dx*dy + j*dx) += srcAddr[i*dy + j] ;
               }
            }
            srcAddr += opCount ;
         }
         ++pmsg ;
      }
      if (colMax) {
         srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi] ;
            for (Index_t i=0; i<dz; ++i) {
               for (Index_t j=0; j<dy; ++j) {
                  (domain.*dest)(dx - 1 + i*dx*dy + j*dx) += srcAddr[i*dy + j] ;
               }
            }
            srcAddr += opCount ;
         }
         ++pmsg ;
      }
   }

   if (rowMin & colMin) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dz; ++i) (domain.*dest)(i*dx*dy) += srcAddr[i] ;
         srcAddr += dz ;
      }
      ++emsg ;
   }

   if (rowMin & planeMin) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dx; ++i) (domain.*dest)(i) += srcAddr[i] ;
         srcAddr += dx ;
      }
      ++emsg ;
   }

   if (colMin & planeMin) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dy; ++i) (domain.*dest)(i*dx) += srcAddr[i] ;
         srcAddr += dy ;
      }
      ++emsg ;
   }

   if (rowMax & colMax) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dz; ++i) (domain.*dest)(dx*dy - 1 + i*dx*dy) += srcAddr[i] ;
         srcAddr += dz ;
      }
      ++emsg ;
   }

   if (rowMax & planeMax) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dx; ++i) (domain.*dest)(dx*(dy-1) + dx*dy*(dz-1) + i) += srcAddr[i] ;
         srcAddr += dx ;
      }
      ++emsg ;
   }

   if (colMax & planeMax) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dy; ++i) (domain.*dest)(dx*dy*(dz-1) + dx - 1 + i*dx) += srcAddr[i] ;
         srcAddr += dy ;
      }
      ++emsg ;
   }

   if (rowMax & colMin) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dz; ++i) (domain.*dest)(dx*(dy-1) + i*dx*dy) += srcAddr[i] ;
         srcAddr += dz ;
      }
      ++emsg ;
   }

   if (rowMin & planeMax) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dx; ++i) (domain.*dest)(dx*dy*(dz-1) + i) += srcAddr[i] ;
         srcAddr += dx ;
      }
      ++emsg ;
   }

   if (colMin & planeMax) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dy; ++i) (domain.*dest)(dx*dy*(dz-1) + i*dx) += srcAddr[i] ;
         srcAddr += dy ;
      }
      ++emsg ;
   }

   if (rowMin & colMax) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dz; ++i) (domain.*dest)(dx - 1 + i*dx*dy) += srcAddr[i] ;
         srcAddr += dz ;
      }
      ++emsg ;
   }

   if (rowMax & planeMin) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dx; ++i) (domain.*dest)(dx*(dy - 1) + i) += srcAddr[i] ;
         srcAddr += dx ;
      }
      ++emsg ;
   }

   if (colMax & planeMin) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dy; ++i) (domain.*dest)(dx - 1 + i*dx) += srcAddr[i] ;
         srcAddr += dy ;
      }
      ++emsg ;
   }

   if (rowMin & colMin & planeMin) {
      Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
      for (Index_t fi=0; fi<xferFields; ++fi) (domain.*fieldData[fi])(0) += comBuf[fi] ;
      ++cmsg ;
   }
   if (rowMin & colMin & planeMax) {
      Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx*dy*(dz - 1) ;
      for (Index_t fi=0; fi<xferFields; ++fi) (domain.*fieldData[fi])(idx) += comBuf[fi] ;
      ++cmsg ;
   }
   if (rowMin & colMax & planeMin) {
      Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx - 1 ;
      for (Index_t fi=0; fi<xferFields; ++fi) (domain.*fieldData[fi])(idx) += comBuf[fi] ;
      ++cmsg ;
   }
   if (rowMin & colMax & planeMax) {
      Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx*dy*(dz - 1) + (dx - 1) ;
      for (Index_t fi=0; fi<xferFields; ++fi) (domain.*fieldData[fi])(idx) += comBuf[fi] ;
      ++cmsg ;
   }
   if (rowMax & colMin & planeMin) {
      Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx*(dy - 1) ;
      for (Index_t fi=0; fi<xferFields; ++fi) (domain.*fieldData[fi])(idx) += comBuf[fi] ;
      ++cmsg ;
   }
   if (rowMax & colMin & planeMax) {
      Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx*dy*(dz - 1) + dx*(dy - 1) ;
      for (Index_t fi=0; fi<xferFields; ++fi) (domain.*fieldData[fi])(idx) += comBuf[fi] ;
      ++cmsg ;
   }
   if (rowMax & colMax & planeMin) {
      Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx*dy - 1 ;
      for (Index_t fi=0; fi<xferFields; ++fi) (domain.*fieldData[fi])(idx) += comBuf[fi] ;
      ++cmsg ;
   }
   if (rowMax & colMax & planeMax) {
      Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx*dy*dz - 1 ;
      for (Index_t fi=0; fi<xferFields; ++fi) (domain.*fieldData[fi])(idx) += comBuf[fi] ;
      ++cmsg ;
   }
}

/******************************************/

void CommSyncPosVel(Domain& domain) {

   if (domain.numRanks() == 1)
      return ;

   // mdmp_wait(-1); removed - LLVM pass intercepts memory read

   bool doRecv = false ; 
   Index_t xferFields = 6 ; /* x, y, z, xd, yd, zd */
   Domain_member fieldData[6] ;
   Index_t maxPlaneComm = xferFields * domain.maxPlaneSize() ;
   Index_t maxEdgeComm  = xferFields * domain.maxEdgeSize() ;
   Index_t pmsg = 0 ; /* plane comm msg */
   Index_t emsg = 0 ; /* edge comm msg */
   Index_t cmsg = 0 ; /* corner comm msg */
   Index_t dx = domain.sizeX() + 1 ;
   Index_t dy = domain.sizeY() + 1 ;
   Index_t dz = domain.sizeZ() + 1 ;
   Real_t *srcAddr ;
   bool rowMin, rowMax, colMin, colMax, planeMin, planeMax ;

   rowMin = rowMax = colMin = colMax = planeMin = planeMax = true ;
   if (domain.rowLoc() == 0) rowMin = false ;
   if (domain.rowLoc() == (domain.tp()-1)) rowMax = false ;
   if (domain.colLoc() == 0) colMin = false ;
   if (domain.colLoc() == (domain.tp()-1)) colMax = false ;
   if (domain.planeLoc() == 0) planeMin = false ;
   if (domain.planeLoc() == (domain.tp()-1)) planeMax = false ;

   fieldData[0] = &Domain::x ;
   fieldData[1] = &Domain::y ;
   fieldData[2] = &Domain::z ;
   fieldData[3] = &Domain::xd ;
   fieldData[4] = &Domain::yd ;
   fieldData[5] = &Domain::zd ;

   if (planeMin | planeMax) {
      Index_t opCount = dx * dy ;

      if (planeMin && doRecv) {
         srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi] ;
            for (Index_t i=0; i<opCount; ++i) (domain.*dest)(i) = srcAddr[i] ;
            srcAddr += opCount ;
         }
         ++pmsg ;
      }
      if (planeMax) {
         srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi] ;
            for (Index_t i=0; i<opCount; ++i) (domain.*dest)(dx*dy*(dz - 1) + i) = srcAddr[i] ;
            srcAddr += opCount ;
         }
         ++pmsg ;
      }
   }

   if (rowMin | rowMax) {
      Index_t opCount = dx * dz ;

      if (rowMin && doRecv) {
         srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi] ;
            for (Index_t i=0; i<dz; ++i) {
               for (Index_t j=0; j<dx; ++j) (domain.*dest)(i*dx*dy + j) = srcAddr[i*dx + j] ;
            }
            srcAddr += opCount ;
         }
         ++pmsg ;
      }
      if (rowMax) {
         srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi] ;
            for (Index_t i=0; i<dz; ++i) {
               for (Index_t j=0; j<dx; ++j) (domain.*dest)(dx*(dy - 1) + i*dx*dy + j) = srcAddr[i*dx + j] ;
            }
            srcAddr += opCount ;
         }
         ++pmsg ;
      }
   }

   if (colMin | colMax) {
      Index_t opCount = dy * dz ;

      if (colMin && doRecv) {
         srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi] ;
            for (Index_t i=0; i<dz; ++i) {
               for (Index_t j=0; j<dy; ++j) (domain.*dest)(i*dx*dy + j*dx) = srcAddr[i*dy + j] ;
            }
            srcAddr += opCount ;
         }
         ++pmsg ;
      }
      if (colMax) {
         srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi] ;
            for (Index_t i=0; i<dz; ++i) {
               for (Index_t j=0; j<dy; ++j) (domain.*dest)(dx - 1 + i*dx*dy + j*dx) = srcAddr[i*dy + j] ;
            }
            srcAddr += opCount ;
         }
         ++pmsg ;
      }
   }

   if (rowMin && colMin && doRecv) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dz; ++i) (domain.*dest)(i*dx*dy) = srcAddr[i] ;
         srcAddr += dz ;
      }
      ++emsg ;
   }

   if (rowMin && planeMin && doRecv) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dx; ++i) (domain.*dest)(i) = srcAddr[i] ;
         srcAddr += dx ;
      }
      ++emsg ;
   }

   if (colMin && planeMin && doRecv) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dy; ++i) (domain.*dest)(i*dx) = srcAddr[i] ;
         srcAddr += dy ;
      }
      ++emsg ;
   }

   if (rowMax && colMax) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dz; ++i) (domain.*dest)(dx*dy - 1 + i*dx*dy) = srcAddr[i] ;
         srcAddr += dz ;
      }
      ++emsg ;
   }

   if (rowMax && planeMax) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dx; ++i) (domain.*dest)(dx*(dy-1) + dx*dy*(dz-1) + i) = srcAddr[i] ;
         srcAddr += dx ;
      }
      ++emsg ;
   }

   if (colMax && planeMax) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dy; ++i) (domain.*dest)(dx*dy*(dz-1) + dx - 1 + i*dx) = srcAddr[i] ;
         srcAddr += dy ;
      }
      ++emsg ;
   }

   if (rowMax && colMin) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dz; ++i) (domain.*dest)(dx*(dy-1) + i*dx*dy) = srcAddr[i] ;
         srcAddr += dz ;
      }
      ++emsg ;
   }

   if (rowMin && planeMax) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dx; ++i) (domain.*dest)(dx*dy*(dz-1) + i) = srcAddr[i] ;
         srcAddr += dx ;
      }
      ++emsg ;
   }

   if (colMin && planeMax) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dy; ++i) (domain.*dest)(dx*dy*(dz-1) + i*dx) = srcAddr[i] ;
         srcAddr += dy ;
      }
      ++emsg ;
   }

   if (rowMin && colMax && doRecv) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dz; ++i) (domain.*dest)(dx - 1 + i*dx*dy) = srcAddr[i] ;
         srcAddr += dz ;
      }
      ++emsg ;
   }

   if (rowMax && planeMin && doRecv) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dx; ++i) (domain.*dest)(dx*(dy - 1) + i) = srcAddr[i] ;
         srcAddr += dx ;
      }
      ++emsg ;
   }

   if (colMax && planeMin && doRecv) {
      srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Domain_member dest = fieldData[fi] ;
         for (Index_t i=0; i<dy; ++i) (domain.*dest)(dx - 1 + i*dx) = srcAddr[i] ;
         srcAddr += dy ;
      }
      ++emsg ;
   }


   if (rowMin && colMin && planeMin && doRecv) {
      Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
      for (Index_t fi=0; fi<xferFields; ++fi) (domain.*fieldData[fi])(0) = comBuf[fi] ;
      ++cmsg ;
   }
   if (rowMin && colMin && planeMax) {
      Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx*dy*(dz - 1) ;
      for (Index_t fi=0; fi<xferFields; ++fi) (domain.*fieldData[fi])(idx) = comBuf[fi] ;
      ++cmsg ;
   }
   if (rowMin && colMax && planeMin && doRecv) {
      Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx - 1 ;
      for (Index_t fi=0; fi<xferFields; ++fi) (domain.*fieldData[fi])(idx) = comBuf[fi] ;
      ++cmsg ;
   }
   if (rowMin && colMax && planeMax) {
      Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx*dy*(dz - 1) + (dx - 1) ;
      for (Index_t fi=0; fi<xferFields; ++fi) (domain.*fieldData[fi])(idx) = comBuf[fi] ;
      ++cmsg ;
   }
   if (rowMax && colMin && planeMin && doRecv) {
      Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx*(dy - 1) ;
      for (Index_t fi=0; fi<xferFields; ++fi) (domain.*fieldData[fi])(idx) = comBuf[fi] ;
      ++cmsg ;
   }
   if (rowMax && colMin && planeMax) {
      Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx*dy*(dz - 1) + dx*(dy - 1) ;
      for (Index_t fi=0; fi<xferFields; ++fi) (domain.*fieldData[fi])(idx) = comBuf[fi] ;
      ++cmsg ;
   }
   if (rowMax && colMax && planeMin && doRecv) {
      Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx*dy - 1 ;
      for (Index_t fi=0; fi<xferFields; ++fi) (domain.*fieldData[fi])(idx) = comBuf[fi] ;
      ++cmsg ;
   }
   if (rowMax && colMax && planeMax) {
      Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx*dy*dz - 1 ;
      for (Index_t fi=0; fi<xferFields; ++fi) (domain.*fieldData[fi])(idx) = comBuf[fi] ;
      ++cmsg ;
   }
}

/******************************************/

void CommMonoQ(Domain& domain)
{
   if (domain.numRanks() == 1)
      return ;

   // mdmp_wait(-1); removed - LLVM pass intercepts memory read

   Index_t xferFields = 3 ; /* delv_xi, delv_eta, delv_zeta */
   Domain_member fieldData[3] ;
   Index_t fieldOffset[3] ;
   Index_t maxPlaneComm = xferFields * domain.maxPlaneSize() ;
   Index_t pmsg = 0 ; /* plane comm msg */
   Index_t dx = domain.sizeX() ;
   Index_t dy = domain.sizeY() ;
   Index_t dz = domain.sizeZ() ;
   Real_t *srcAddr ;
   bool rowMin, rowMax, colMin, colMax, planeMin, planeMax ;

   rowMin = rowMax = colMin = colMax = planeMin = planeMax = true ;
   if (domain.rowLoc() == 0) rowMin = false ;
   if (domain.rowLoc() == (domain.tp()-1)) rowMax = false ;
   if (domain.colLoc() == 0) colMin = false ;
   if (domain.colLoc() == (domain.tp()-1)) colMax = false ;
   if (domain.planeLoc() == 0) planeMin = false ;
   if (domain.planeLoc() == (domain.tp()-1)) planeMax = false ;

   fieldData[0] = &Domain::delv_xi ;
   fieldData[1] = &Domain::delv_eta ;
   fieldData[2] = &Domain::delv_zeta ;
   fieldOffset[0] = domain.numElem() ;
   fieldOffset[1] = domain.numElem() ;
   fieldOffset[2] = domain.numElem() ;

   if (planeMin | planeMax) {
      Index_t opCount = dx * dy ;

      if (planeMin) {
         srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi] ;
            for (Index_t i=0; i<opCount; ++i) (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i] ;
            srcAddr += opCount ;
            fieldOffset[fi] += opCount ;
         }
         ++pmsg ;
      }
      if (planeMax) {
         srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi] ;
            for (Index_t i=0; i<opCount; ++i) (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i] ;
            srcAddr += opCount ;
            fieldOffset[fi] += opCount ;
         }
         ++pmsg ;
      }
   }

   if (rowMin | rowMax) {
      Index_t opCount = dx * dz ;

      if (rowMin) {
         srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi] ;
            for (Index_t i=0; i<opCount; ++i) (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i] ;
            srcAddr += opCount ;
            fieldOffset[fi] += opCount ;
         }
         ++pmsg ;
      }
      if (rowMax) {
         srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi] ;
            for (Index_t i=0; i<opCount; ++i) (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i] ;
            srcAddr += opCount ;
            fieldOffset[fi] += opCount ;
         }
         ++pmsg ;
      }
   }
   if (colMin | colMax) {
      Index_t opCount = dy * dz ;

      if (colMin) {
         srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi] ;
            for (Index_t i=0; i<opCount; ++i) (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i] ;
            srcAddr += opCount ;
            fieldOffset[fi] += opCount ;
         }
         ++pmsg ;
      }
      if (colMax) {
         srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi] ;
            for (Index_t i=0; i<opCount; ++i) (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i] ;
            srcAddr += opCount ;
            fieldOffset[fi] += opCount ;
         }
         ++pmsg ;
      }
   }
}

#endif
