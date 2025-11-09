#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <cmath>
#include <float.h>
#include <mpi.h>
#include <pnetcdf.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

//#define PRINT_STATUS
#ifdef PRINT_STATUS
#define printStatus() printf("Status %d: %d\n", rank, statcnt++ )
#else
#define printStatus()
#endif


void checkError( int rank, int status, int line ) {
  if ( status != NC_NOERR ) {
    printf("error from rank %d at line %d: %s\n", rank, line, ncmpi_strerror(status) );
    exit(1);
  }
}

int main(int argc, char **argv){
  int rank, ncpus;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &ncpus );
  int noCells;
  long noCellsLong;
  int ncId0, ncId1, dimId, varId, status;
  MPI_Offset size, start, count;
  if ( argc != 3 ) {
    if ( rank == 0 ) {
      if ( argc < 3 ) printf( "\nNot enough arguments. Usage: >> 'convert_grid_file FILENAME.Netcdf NEWFILE.Netcdf' <<\n\n" );
      if ( argc > 3 ) printf( "\nToo many arguments. Usage: >> 'convert_grid_file FILENAME.Netcdf NEWFILE.Netcdf' <<\n\n" );
    }
    MPI_Finalize();
    return 0;
  }
  const size_t flen = strlen(argv[1]);
  const size_t flen2 = strlen(argv[2]);
  char* fileName = malloc(flen+1);
  char* newFileName = malloc(flen2+1);
  strcpy( fileName, argv[1] );
  strcpy( newFileName, argv[2] );

  if ( access( fileName, F_OK) == -1 ) {
    printf( "\nFile not found: %s \n\n", fileName );
    MPI_Finalize();
    return 1;
  }

  if( access( newFileName, F_OK ) != -1 ) {
    printf( "\nError: file already exists: %s \n\n", newFileName );
    MPI_Finalize();
    return 1;
  }

  char* exts = strrchr(fileName, '.' );
  if ( !exts || strcmp( exts, ".Netcdf" ) != 0 ) {
    printf( "\nInvalid extension for file: %s \n\n", fileName );
    MPI_Finalize();
    return 1;
  }
  exts = strrchr(newFileName, '.' );
  if ( !exts || strcmp( exts, ".Netcdf" ) != 0 ) {
    printf( "\nInvalid extension for file: %s \n\n", newFileName );
    MPI_Finalize();
    return 1;
  }

  MPI_Barrier(MPI_COMM_WORLD);


#ifdef PRINT_STATUS
  int statcnt = 0;
#endif
  int nDim = -1;
  int maxChilds = -1;
  int maxNghbrs = -1;
  int noPartitionCells = -1;
  long noPartitionCellsLong = -1;
  long noLeafCells = -1;
  int noLocalCells = -1;
  int localCellOffset = -1;
  int minLevel = -1;
  int maxLevel = -1;
  int maxUniformRefinementLevel = -1;
  double totalWorkload = -1.;
  double lengthLevel0 = -1.;
  double maxCpus = -1.;
  int minMaxLevelInFile = 1;
  int lengthLevel0InFile = 1;
  int centerOfGravityInFile = 1;
  double centerOfGravity[3] = {.0,.0,.0};
  double bbox[6] = {.0,.0,.0,.0,.0,.0};
  int decisiveDirection = -1;
  double reductionFactor = -1.0;

  int noLocalMinLevelCells = 0;
  int minLevelCellOffset = -1;
  int noMinLevelCells = -1;
  long noMinLevelCellsLong = -1;

  long noPartitionLevelAncestors = 0;
  int maxPartitionLevelShift = 0;
  int globalTimeStep = 0;
  
  int* offsets = (int*)malloc((ncpus+1)*sizeof(int));
  int* domainOffsets = (int*)malloc((ncpus+1)*sizeof(int));
  int* partitionCellsGlobalId = NULL;
  long* partitionCellsGlobalIdLong = NULL;
  double* partitionCellsWorkload = NULL;
  int* partitionCellsLevelDiff = NULL;

  int* parentId = NULL;
  int* noChildIds = NULL;
  int* childIds = NULL;
  int* position = NULL;
  int* nghbrIds = NULL;
  int* level = NULL;
  double* coordinates = NULL;
  
  unsigned char* cellInfo = NULL;
  long* minLevelTreeId = NULL;
  long* minLevelNghbrIds = NULL;
  double* minLevelCoords = NULL;

  MPI_Barrier(MPI_COMM_WORLD);
  const double time0 = MPI_Wtime();
  int status0 = ncmpi_open(MPI_COMM_WORLD, fileName, NC_NOWRITE|NC_64BIT_DATA, MPI_INFO_NULL, &ncId0);
  if (status0 == NC_NOERR) {
    if ( rank == 0 ) {
      printf("Loading file %s using %d cores\n", fileName, ncpus );
    }
    
    printStatus();

    start = 0;
    count = 1;
    varId = -1;
    status = ncmpi_inq_varid(ncId0, "noCells", &varId);
    checkError(rank,status,__LINE__);
    status = ncmpi_get_vara_int_all(ncId0, varId, &start, &count, &noCells);
    checkError(rank,status,__LINE__);
    noCellsLong = (long)noCells;
    
    varId = -1;
    status = ncmpi_inq_varid(ncId0, "minLevel", &varId);
    if ( varId == -1 ) {
      minMaxLevelInFile = 0;
    } else {
      checkError(rank,status,__LINE__);
      status = ncmpi_get_vara_int_all(ncId0, varId, &start, &count, &minLevel);
      checkError(rank,status,__LINE__);
    }

    varId = -1;
    status = ncmpi_inq_varid(ncId0, "maxLevel", &varId);
    if ( varId == -1 ) {
      minMaxLevelInFile = 0;
    } else {
      checkError(rank,status,__LINE__);
      status = ncmpi_get_vara_int_all(ncId0, varId, &start, &count, &maxLevel);
      checkError(rank,status,__LINE__);
    }

    varId = -1;
    status = ncmpi_inq_varid(ncId0, "lengthLevel0", &varId);
    if ( varId == -1 ) {
      lengthLevel0InFile = 0;
    } else {
      checkError(rank,status,__LINE__);
      status = ncmpi_get_vara_double_all(ncId0, varId, &start, &count, &lengthLevel0);
      checkError(rank,status,__LINE__);
    }

    varId = -1;
    status = ncmpi_inq_varid(ncId0, "maxCpus", &varId);
    if ( varId != -1 ) {
      checkError(rank,status,__LINE__);
      status = ncmpi_get_vara_double_all(ncId0, varId, &start, &count, &maxCpus);
      checkError(rank,status,__LINE__);
    }

    varId = -1;
    status = ncmpi_inq_varid(ncId0, "level_2", &varId);
    if ( varId == -1 ) {
      nDim = 2;
    } else {
      nDim = 3;
    }
    maxChilds = (int)pow(2,nDim);
    maxNghbrs = 2*nDim;

    varId = -1;
    status = ncmpi_inq_varid(ncId0, "centerOfGravity", &varId);
    count = nDim;
    if ( varId == -1 ) {
      centerOfGravityInFile = 0;
    } else {
      checkError(rank,status,__LINE__);
      status = ncmpi_get_vara_double_all(ncId0, varId, &start, &count, &centerOfGravity[0]);
      checkError(rank,status,__LINE__);
    }

    printStatus();
    
    varId = -1;
    status = ncmpi_inq_varid(ncId0, "minCellsId", &varId);
    checkError(rank,status,__LINE__);
    status = ncmpi_inq_vardimid(ncId0, varId, &dimId);
    checkError(rank,status,__LINE__);
    status = ncmpi_inq_dimlen(ncId0, dimId, &size);
    checkError(rank,status,__LINE__);
    noPartitionCells = (int)size;
    noPartitionCellsLong = (long)size;
    
    printStatus();

    if ( rank == 0 ) {
      partitionCellsGlobalId = (int*)malloc(noPartitionCells*sizeof(int));
      partitionCellsGlobalIdLong = (long*)malloc(noPartitionCells*sizeof(long));
      partitionCellsWorkload = (double*)malloc(noPartitionCells*sizeof(double));
      partitionCellsLevelDiff = (int*)malloc(noPartitionCells*sizeof(int));
      start = 0;
      count = noPartitionCells;
      status = ncmpi_inq_varid(ncId0, "minCellsId", &varId);
      checkError(rank,status,__LINE__);
      status = ncmpi_get_vara_int_all(ncId0, varId, &start, &count, partitionCellsGlobalId);
      checkError(rank,status,__LINE__);
      status = ncmpi_inq_varid(ncId0, "minCellsWorkLoad", &varId);
      checkError(rank,status,__LINE__);
      status = ncmpi_get_vara_double_all(ncId0, varId, &start, &count, partitionCellsWorkload);
      checkError(rank,status,__LINE__);
      status = ncmpi_inq_varid(ncId0, "minCellsLvlDiff", &varId);
      checkError(rank,status,__LINE__);
      status = ncmpi_get_vara_int_all(ncId0, varId, &start, &count, partitionCellsLevelDiff);
      checkError(rank,status,__LINE__);

      maxPartitionLevelShift = 0;
      for ( int i = 0; i < noPartitionCells; i++ ) {
        maxPartitionLevelShift = MAX( maxPartitionLevelShift, partitionCellsLevelDiff[i] );
      }

      for ( int i = 0; i < noPartitionCells; i++ ) {
        partitionCellsGlobalIdLong[i] = (long)partitionCellsGlobalId[i];
      }
      
      free(partitionCellsGlobalId); partitionCellsGlobalId = NULL;
      free(partitionCellsLevelDiff); partitionCellsLevelDiff = NULL;

      totalWorkload = 0.0;
      for ( int i = 0; i < noPartitionCells; i++ ) {
        totalWorkload += partitionCellsWorkload[i];
      }
      const double avgWL = totalWorkload / ((double)ncpus);
      offsets[0] = 0;
      double counter = 0.0;
      for ( int i = 0; i < noPartitionCells; i++ ) {
        int f0 = (int)(counter/avgWL);
        int f1 = (int)((counter+partitionCellsWorkload[i])/avgWL);
        if ( f1 > f0 && f1 < ncpus ) {
          offsets[f1] = i+1;
        }
        counter += partitionCellsWorkload[i];
      }
      offsets[ncpus] = noPartitionCells;
      domainOffsets[0] = 0;
      for ( int i = 1; i < ncpus; i++ ) {
        domainOffsets[i] = (int)partitionCellsGlobalIdLong[ offsets[i] ];
      }
      domainOffsets[ncpus] = noCells;
    }
    else {
      double fltDummy;
      int intDummy;
      start = 0;
      count = 0;
      status = ncmpi_inq_varid(ncId0, "minCellsId", &varId);
      checkError(rank,status,__LINE__);
      status = ncmpi_get_vara_int_all(ncId0, varId, &start, &count, &intDummy);
      checkError(rank,status,__LINE__);
      status = ncmpi_inq_varid(ncId0, "minCellsWorkLoad", &varId);
      checkError(rank,status,__LINE__);
      status = ncmpi_get_vara_double_all(ncId0, varId, &start, &count, &fltDummy);
      checkError(rank,status,__LINE__);
      status = ncmpi_inq_varid(ncId0, "minCellsLvlDiff", &varId);
      checkError(rank,status,__LINE__);
      status = ncmpi_get_vara_int_all(ncId0, varId, &start, &count, &intDummy);
      checkError(rank,status,__LINE__);
    }

    printStatus();
    
    MPI_Bcast(offsets, ncpus+1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(domainOffsets, ncpus+1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&totalWorkload, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxPartitionLevelShift, 1, MPI_INT, 0, MPI_COMM_WORLD);

    printStatus();

    noLocalCells = domainOffsets[rank+1] - domainOffsets[rank];
    localCellOffset = 0;
    MPI_Exscan( &noLocalCells, &localCellOffset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    
    printStatus();

    start = localCellOffset;
    count = noLocalCells;
    status = ncmpi_inq_varid(ncId0, "parentId", &varId);
    checkError(rank,status,__LINE__);

    printStatus();

    parentId = (int*)malloc(noLocalCells*sizeof(int));
    status = ncmpi_get_vara_int_all(ncId0, varId, &start, &count, parentId);
    checkError(rank,status,__LINE__);

    printStatus();

    int noLocalCellsMissing = 0;
    int maxGlobalCellsMissing = 0;
    for ( int i = 0; i < noLocalCells; i++ ) {
      if ( parentId[i] > -1 && (parentId[i] < localCellOffset || parentId[i] >= localCellOffset+noLocalCells) ) {
        noLocalCellsMissing++;
      }
    }
    MPI_Allreduce(&noLocalCellsMissing, &maxGlobalCellsMissing, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    int* missingParents = (int*)malloc(noLocalCellsMissing*sizeof(int));
    int pcnt = 0;
    for ( int i = 0; i < noLocalCells; i++ ) {
      if ( parentId[i] > -1 && (parentId[i] < localCellOffset || parentId[i] >= localCellOffset+noLocalCells) ) {
        missingParents[ pcnt++ ] = parentId[i];
      }
    }

    childIds = (int*)malloc((noLocalCells+noLocalCellsMissing)*sizeof(int));
    position = (int*)malloc(noLocalCells*sizeof(int));
    for ( int i = 0; i < noLocalCells; i++ ) {
      position[i] = 0;
    }
    for ( int k = 0; k < maxChilds; k++ ) {
      char str[11];
      sprintf(str, "childIds_%d", k);
      status = ncmpi_inq_varid(ncId0, str, &varId);
      checkError(rank,status,__LINE__);
      start = localCellOffset;
      count = noLocalCells;
      status = ncmpi_get_vara_int_all(ncId0, varId, &start, &count, childIds);
      checkError(rank,status,__LINE__);
      pcnt = 0;
      for ( int k = 0; k < noLocalCellsMissing; k++ ) {
        start = missingParents[ pcnt ];
        count = 1;
        status = ncmpi_get_vara_int_all(ncId0, varId, &start, &count, &childIds[noLocalCells+pcnt]);
        pcnt++;
        checkError(rank,status,__LINE__);
      }
      for ( int k = 0; k < maxGlobalCellsMissing-noLocalCellsMissing; k++ ) {
        int intDummy;
        start = 0;
        count = 0;
        status = ncmpi_get_vara_int_all(ncId0, varId, &start, &count, &intDummy);
        checkError(rank,status,__LINE__);
      }
      pcnt = 0;
      for ( int i = 0; i < noLocalCells; i++ ) {
        if ( parentId[i] > -1 ) {
          int child = -1;
          if ( (parentId[i] >= localCellOffset && parentId[i] < localCellOffset+noLocalCells) ) {
            int pid = parentId[i] - localCellOffset;
            child = childIds[pid];
          } else {
            child = childIds[noLocalCells+pcnt];
            pcnt++;
          }
          if ( child == localCellOffset+i ) { 
            position[i] = k;
          }
        }
      }
    }
    free( childIds ); childIds = NULL;

    printStatus();

    start = localCellOffset;
    count = noLocalCells;

    noChildIds = (int*)malloc(noLocalCells*sizeof(int));
    status = ncmpi_inq_varid(ncId0, "noChildIds", &varId);
    checkError(rank,status,__LINE__);
    status = ncmpi_get_vara_int_all(ncId0, varId, &start, &count, noChildIds);
    checkError(rank,status,__LINE__);

    printStatus();

    cellInfo = (unsigned char*)malloc(noLocalCells*sizeof(unsigned char));

    noLocalMinLevelCells = 0;
    noLeafCells = 0;
    for ( int i = 0; i < noLocalCells; i++ ) {
      unsigned int noChilds = (unsigned int)noChildIds[i];
      unsigned int isMinLvl = ( parentId[i] < 0 );
      unsigned int pos = (unsigned int)position[i];
      unsigned int tmpBit = noChilds | (pos << 4) | (isMinLvl << 7);
      cellInfo[i] = (unsigned char)tmpBit;
      if ( parentId[i] < 0 ) {
        noLocalMinLevelCells++;
      }
      if ( noChildIds[i] == 0 ) {
        noLeafCells++;
      }
    }
    

    free( noChildIds ); noChildIds = NULL;
    free( position ); position = NULL;
    
    printStatus();

        
    minLevelTreeId = (long*)malloc(noLocalMinLevelCells*sizeof(long));
    minLevelNghbrIds = (long*)malloc(maxNghbrs*noLocalMinLevelCells*sizeof(long));
    minLevelCoords = (double*)malloc(nDim*noLocalMinLevelCells*sizeof(double));
    
    minLevelCellOffset = 0;
    noMinLevelCells = 0;
    MPI_Exscan( &noLocalMinLevelCells, &minLevelCellOffset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    MPI_Allreduce(&noLocalMinLevelCells, &noMinLevelCells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &noLeafCells, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    noMinLevelCellsLong = (long)noMinLevelCells;

    printStatus();

    level = (int*)malloc(noLocalCells*sizeof(int));
    coordinates = (double*)malloc(noLocalCells*sizeof(double));
    nghbrIds = (int*)malloc(noLocalCells*sizeof(int));

    status = ncmpi_inq_varid(ncId0, "level_0", &varId);
    checkError(rank,status,__LINE__);
    status = ncmpi_get_vara_int_all(ncId0, varId, &start, &count, level);
    checkError(rank,status,__LINE__);

    printStatus();

    if ( !minMaxLevelInFile ) {
      minLevel = 9999;
      maxLevel = -1;
      for ( int i = 0; i < noLocalCells; i++ ) {
        minLevel = MIN( minLevel, level[i] );
        maxLevel = MAX( maxLevel, level[i] );
      }
      MPI_Allreduce(MPI_IN_PLACE, &minLevel, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &maxLevel, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }

    maxUniformRefinementLevel = maxLevel;
    for ( int i = 0; i < noLocalCells; i++ ) {
      unsigned int cinfo = (unsigned int)cellInfo[i];
      unsigned int noChilds = cinfo & 15u;
      if ( noChilds == 0 ) {
        maxUniformRefinementLevel = MIN( maxUniformRefinementLevel, level[i] );
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &maxUniformRefinementLevel, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    for ( int k = 0; k < nDim; k++ ) {
      bbox[k] = DBL_MAX;
      bbox[nDim+k] = -DBL_MAX;
    }

    if ( !lengthLevel0InFile ) {
      status = ncmpi_inq_varid(ncId0, "coordinates_0", &varId);
      checkError(rank,status,__LINE__);
      status = ncmpi_get_vara_double_all(ncId0, varId, &start, &count, coordinates);
      checkError(rank,status,__LINE__);
      status = ncmpi_inq_varid(ncId0, "nghbrIds_0", &varId);
      checkError(rank,status,__LINE__);
      status = ncmpi_get_vara_int_all(ncId0, varId, &start, &count, nghbrIds);
      checkError(rank,status,__LINE__);
      int nghbrsFound = 0;
      int id = 0;
      while ( !nghbrsFound && id < noLocalCells ) {
        int nghbr = nghbrIds[id];
        if ( nghbr > -1 && nghbr >= localCellOffset && nghbr < localCellOffset+noLocalCells ) {
          nghbr -= localCellOffset;
          if ( level[id] == level[nghbr] ) {
            double x1 = coordinates[id];
            double x0 = coordinates[nghbr];
            lengthLevel0 = (x1-x0)*pow(2.0,(double)level[id]);
            nghbrsFound = 1;
          }
        } else {
          id++;
        }
      }
      if ( !nghbrsFound ) {
        printf( "This was not supposed to happen.." );
        MPI_Finalize();
        return 1;
      }
      MPI_Allreduce(MPI_IN_PLACE, &lengthLevel0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      lengthLevel0 /= ((double)ncpus);
    }

    for ( int k = 0; k < nDim; k++ ) {
      char str[14];
      sprintf(str, "coordinates_%d", k);
      status = ncmpi_inq_varid(ncId0, str, &varId);
      checkError(rank,status,__LINE__);
      status = ncmpi_get_vara_double_all(ncId0, varId, &start, &count, coordinates);
      checkError(rank,status,__LINE__);
      int cnt = 0;
      for ( int i = 0; i < noLocalCells; i++ ) {
        if ( parentId[i] > -1 ) continue;
        double dx = 0.5 * lengthLevel0 / pow(2.0,(double)level[i]) ;
        bbox[k] = MIN( bbox[k], (coordinates[i] - dx) );
        bbox[nDim+k] = MAX( bbox[nDim+k], coordinates[i] + dx );
        if ( parentId[i] < 0 ) {
          minLevelCoords[cnt*nDim+k] = coordinates[i];
          cnt++;
        }
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &bbox[0], nDim, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &bbox[nDim], nDim, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if ( !centerOfGravityInFile ) {
      for ( int k = 0; k < nDim; k++ ) {
        centerOfGravity[k] = 0.5*( bbox[k] + bbox[nDim+k] );
        long ncells = (long)( ( bbox[nDim+k] - bbox[k] ) * pow(2.0,(double)minLevel) / lengthLevel0);
        if ( ncells % 2l == 1 ) {
          centerOfGravity[k] += 0.5 * lengthLevel0 / pow(2.0,(double)minLevel);
        }
      }
    }

    double maxExt = 0.0;
    for ( int k = 0; k < nDim; k++ ) {
      if ( bbox[nDim+k] - bbox[k] > maxExt ) {
        decisiveDirection = k;
        maxExt = bbox[nDim+k] - bbox[k];
      }
    }


    printStatus();


    for ( int k = 0; k < maxNghbrs; k++ ) {
      char str[11];
      sprintf(str, "nghbrIds_%d", k);
      status = ncmpi_inq_varid(ncId0, str, &varId);
      checkError(rank,status,__LINE__);
      status = ncmpi_get_vara_int_all(ncId0, varId, &start, &count, nghbrIds);
      checkError(rank,status,__LINE__);
      
      int cnt = 0;
      for ( int i = 0; i < noLocalCells; i++ ) {
        if ( parentId[i] < 0 ) {
          minLevelNghbrIds[cnt*maxNghbrs+k] = nghbrIds[i];
          cnt++;
        }
      }
    }

    free(level); level = NULL;
    free(coordinates); coordinates = NULL;
    free(nghbrIds); nghbrIds = NULL;

    printStatus();

    int cnt = 0;
    for ( int i = 0; i < noLocalCells; i++ ) {
      if ( parentId[i] < 0 ) {
        long treeId = 0;
        long bitCount = 0;
        long coord[3];
        for ( int k = 0; k < nDim; k++ ) {
          coord[k] = (long)( ( minLevelCoords[cnt*nDim+k] - centerOfGravity[k] + 0.5 * lengthLevel0 ) * pow(2.0,(double)minLevel) / lengthLevel0 );
        }
        for ( long lvl = minLevel-1; lvl >= 0; lvl-- ) {
          for ( long k = 0; k < nDim; k++ ) {
            long bit = ( coord[k] / ((long)pow(2.0,lvl)) ) % 2;
            treeId |= bit << bitCount; 
            bitCount++;
          }
        }
        minLevelTreeId[cnt] = treeId;
        cnt++;
      }
    }

    free(minLevelCoords); minLevelCoords = NULL;
    free( parentId ); parentId = NULL;

    printStatus();


    int noLocalPartitionCells = offsets[rank+1] - offsets[rank];
    int localPartitionCellOffset = offsets[rank];
    partitionCellsGlobalId = (int*)malloc(noLocalPartitionCells*sizeof(int));
    start = localPartitionCellOffset;
    count = noLocalPartitionCells;
    status = ncmpi_inq_varid(ncId0, "minCellsId", &varId);
    checkError(rank,status,__LINE__);
    status = ncmpi_get_vara_int_all(ncId0, varId, &start, &count, partitionCellsGlobalId);
    checkError(rank,status,__LINE__);


    long noPartitionLevelAncestorsLocal = noLocalCells;
    for ( int i = 0; i < noLocalPartitionCells; i++ ) {
      int localId = partitionCellsGlobalId[i] - localCellOffset;
      const unsigned int childCnt = ((unsigned int)cellInfo[localId]) & 15u;
      long noFlagged = 1;
      int childId = localId;
      int totalNoChilds = childCnt;
      while ( totalNoChilds ) {
        childId++;
        if ( childId >= noLocalCells ) {
          printf( "child out of range: %d %d", childId, noLocalCells);
          MPI_Finalize();
          return 1;
        }
        noFlagged++;
        totalNoChilds--;
        totalNoChilds += ((unsigned int)cellInfo[childId]) & 15u;
      }
      noPartitionLevelAncestorsLocal -= noFlagged;
    }
    MPI_Allreduce(&noPartitionLevelAncestorsLocal, &noPartitionLevelAncestors, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    free(partitionCellsGlobalId); partitionCellsGlobalId = NULL;
    printStatus();
  }
  else {
    if ( rank == 0 ) {
      printf( "err" );
      MPI_Finalize();
      return 1;
    }
  }
  status0 = ncmpi_close(ncId0);
  if (status0 != NC_NOERR) {
    printf( "err" );
    MPI_Finalize();
    return 1;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  const double time1 = MPI_Wtime();

  int status1 = ncmpi_create(MPI_COMM_WORLD, newFileName, NC_WRITE|NC_64BIT_DATA|NC_NOCLOBBER, MPI_INFO_NULL, &ncId1);
  if (status1 == NC_NOERR) {

    if ( rank == 0 ) {
      printf("Writing converted file %s using %d cores\n", newFileName, ncpus );
    }
    
    printStatus();
    
    status = ncmpi_def_dim(ncId1, "dim0", noPartitionCells, &dimId);
    checkError(rank,status,__LINE__);
    
    status = ncmpi_def_dim(ncId1, "dim1", noMinLevelCells, &dimId);
    checkError(rank,status,__LINE__);

    status = ncmpi_def_dim(ncId1, "dim2", maxNghbrs*noMinLevelCells, &dimId);
    checkError(rank,status,__LINE__);

    status = ncmpi_def_dim(ncId1, "dim3", noCells, &dimId);
    checkError(rank,status,__LINE__);
    

    status = ncmpi_put_att(ncId1, NC_GLOBAL, "nDim", NC_INT, 1, &nDim);
    checkError(rank,status,__LINE__);

    int noBlocks = 1;
    status = ncmpi_put_att(ncId1, NC_GLOBAL, "noBlocks", NC_INT, 1, &noBlocks);
    checkError(rank,status,__LINE__);

    status = ncmpi_put_att(ncId1, NC_GLOBAL, "globalTimeStep", NC_INT, 1, &globalTimeStep);
    checkError(rank,status,__LINE__);

    status = ncmpi_put_att(ncId1, NC_GLOBAL, "noCells", NC_INT64, 1, &noCellsLong);
    checkError(rank,status,__LINE__);

    status = ncmpi_put_att(ncId1, NC_GLOBAL, "noLeafCells", NC_INT64, 1, &noLeafCells);
    checkError(rank,status,__LINE__);

    status = ncmpi_put_att(ncId1, NC_GLOBAL, "noMinLevelCells", NC_INT64, 1, &noMinLevelCellsLong);
    checkError(rank,status,__LINE__);

    status = ncmpi_put_att(ncId1, NC_GLOBAL, "noPartitionCells", NC_INT64, 1, &noPartitionCellsLong);
    checkError(rank,status,__LINE__);

    status = ncmpi_put_att(ncId1, NC_GLOBAL, "noPartitionLevelAncestors", NC_INT64, 1, &noPartitionLevelAncestors);
    checkError(rank,status,__LINE__);

    status = ncmpi_put_att(ncId1, NC_GLOBAL, "minLevel", NC_INT, 1, &minLevel);
    checkError(rank,status,__LINE__);

    status = ncmpi_put_att(ncId1, NC_GLOBAL, "maxLevel", NC_INT, 1, &maxLevel);
    checkError(rank,status,__LINE__);

    status = ncmpi_put_att(ncId1, NC_GLOBAL, "maxUniformRefinementLevel", NC_INT, 1, &maxUniformRefinementLevel);
    checkError(rank,status,__LINE__);

    status = ncmpi_put_att(ncId1, NC_GLOBAL, "maxPartitionLevelShift", NC_INT, 1, &maxPartitionLevelShift);
    checkError(rank,status,__LINE__);

    status = ncmpi_put_att(ncId1, NC_GLOBAL, "lengthLevel0", NC_DOUBLE, 1, &lengthLevel0);
    checkError(rank,status,__LINE__);

    status = ncmpi_put_att(ncId1, NC_GLOBAL, "centerOfGravity", NC_DOUBLE, nDim, &centerOfGravity);
    checkError(rank,status,__LINE__);

    status = ncmpi_put_att(ncId1, NC_GLOBAL, "boundingBox", NC_DOUBLE, 2*nDim, &bbox);
    checkError(rank,status,__LINE__);

    status = ncmpi_put_att(ncId1, NC_GLOBAL, "reductionFactor", NC_DOUBLE, 1, &reductionFactor);
    checkError(rank,status,__LINE__);

    status = ncmpi_put_att(ncId1, NC_GLOBAL, "decisiveDirection", NC_INT, 1, &decisiveDirection);
    checkError(rank,status,__LINE__);

    status = ncmpi_put_att(ncId1, NC_GLOBAL, "totalWorkload", NC_DOUBLE, 1, &totalWorkload);
    checkError(rank,status,__LINE__);

    long partitionCellOffspringThreshold = -1;
    double partitionCellWorkloadThreshold = -1.;
    status = ncmpi_put_att(ncId1, NC_GLOBAL, "partitionCellOffspringThreshold", NC_INT, 1, &partitionCellOffspringThreshold);
    checkError(rank,status,__LINE__);
    status = ncmpi_put_att(ncId1, NC_GLOBAL, "partitionCellWorkloadThreshold", NC_DOUBLE, 1, &partitionCellWorkloadThreshold);
    checkError(rank,status,__LINE__);

    status = ncmpi_put_att(ncId1, NC_GLOBAL, "maxNoBalancedCPUs", NC_DOUBLE, 1, &maxCpus);
    checkError(rank,status,__LINE__);

    printStatus();
    
    status = ncmpi_inq_dimid(ncId1, "dim0", &dimId);
    checkError(rank,status,__LINE__);
    status = ncmpi_def_var(ncId1, "partitionCellsGlobalId", NC_INT64, 1, &dimId, &varId);
    checkError(rank,status,__LINE__);
    status = ncmpi_def_var(ncId1, "partitionCellsWorkload", NC_DOUBLE, 1, &dimId, &varId);
    checkError(rank,status,__LINE__);
    
    status = ncmpi_inq_dimid(ncId1, "dim1", &dimId);
    checkError(rank,status,__LINE__);
    status = ncmpi_def_var(ncId1, "minLevelCellsTreeId", NC_INT64, 1, &dimId, &varId);
    checkError(rank,status,__LINE__);
    
    status = ncmpi_inq_dimid(ncId1, "dim2", &dimId);
    checkError(rank,status,__LINE__);
    status = ncmpi_def_var(ncId1, "minLevelCellsNghbrIds", NC_INT64, 1, &dimId, &varId);
    checkError(rank,status,__LINE__);

    status = ncmpi_inq_dimid(ncId1, "dim3", &dimId);
    checkError(rank,status,__LINE__);
    status = ncmpi_def_var(ncId1, "cellInfo", NC_UBYTE, 1, &dimId, &varId);
    checkError(rank,status,__LINE__);

    status = ncmpi_enddef( ncId1 );
    checkError(rank,status,__LINE__);

    printStatus();

    if (rank == 0 ) {
      start = 0;
      count = noPartitionCells;
      
      status = ncmpi_inq_varid(ncId1, "partitionCellsGlobalId", &varId);
      checkError(rank,status,__LINE__);
      status = ncmpi_put_vara_long_all(ncId1, varId, &start, &count, partitionCellsGlobalIdLong);
      checkError(rank,status,__LINE__);
      
      status = ncmpi_inq_varid(ncId1, "partitionCellsWorkload", &varId);
      checkError(rank,status,__LINE__);
      status = ncmpi_put_vara_double_all(ncId1, varId, &start, &count, partitionCellsWorkload);
      checkError(rank,status,__LINE__);
    } else {
      double fltDummy;
      long longDummy;
      start = 0;
      count = 0;
      
      status = ncmpi_inq_varid(ncId1, "partitionCellsGlobalId", &varId);
      checkError(rank,status,__LINE__);
      status = ncmpi_put_vara_long_all(ncId1, varId, &start, &count, &longDummy);
      checkError(rank,status,__LINE__);
      
      status = ncmpi_inq_varid(ncId1, "partitionCellsWorkload", &varId);
      checkError(rank,status,__LINE__);
      status = ncmpi_put_vara_double_all(ncId1, varId, &start, &count, &fltDummy);
      checkError(rank,status,__LINE__);
    }

    printStatus();

    //printf("offset %d %d %d \n", minLevelCellOffset, minLevelCellOffset+noLocalMinLevelCells, noMinLevelCells );

    start = minLevelCellOffset;
    count = noLocalMinLevelCells;
      
    status = ncmpi_inq_varid(ncId1, "minLevelCellsTreeId", &varId);
    checkError(rank,status,__LINE__);
    status = ncmpi_put_vara_long_all(ncId1, varId, &start, &count, minLevelTreeId);
    checkError(rank,status,__LINE__);

    printStatus();

    start = maxNghbrs*minLevelCellOffset;
    count = maxNghbrs*noLocalMinLevelCells;
      
    status = ncmpi_inq_varid(ncId1, "minLevelCellsNghbrIds", &varId);
    checkError(rank,status,__LINE__);
    status = ncmpi_put_vara_long_all(ncId1, varId, &start, &count, minLevelNghbrIds);
    checkError(rank,status,__LINE__);

    printStatus();
    
    start = localCellOffset;
    count = noLocalCells;
      
    status = ncmpi_inq_varid(ncId1, "cellInfo", &varId);
    checkError(rank,status,__LINE__);
    status = ncmpi_put_vara_uchar_all(ncId1, varId, &start, &count, cellInfo);
    checkError(rank,status,__LINE__);

    printStatus();

  }
  else {
    if ( status1 == NC_EEXIST ) {
      printf( "Error: file already exists: %s \n", newFileName );
    } else {
      printf( "Error creating file %s \n", newFileName );
    }
    MPI_Finalize();
    return 1;
  }
  status1 = ncmpi_close(ncId1);
  if (status1 != NC_NOERR) {
    printf( "err" );
    MPI_Finalize();
    return 1;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  const double time2 = MPI_Wtime();

  if ( rank == 0 ) {
    printf( "Time for reading / writing in seconds: %f / %f\n", time1-time0, time2-time1 );
  }

  if ( rank == 0 ) {
    free(partitionCellsWorkload); partitionCellsWorkload = NULL;
    free(partitionCellsGlobalIdLong); partitionCellsGlobalIdLong = NULL;
  }
  free(cellInfo); cellInfo = NULL;
  free(minLevelTreeId); minLevelTreeId = NULL;
  free(minLevelNghbrIds); minLevelNghbrIds = NULL;


  free(newFileName); newFileName = NULL;
  free(fileName); fileName = NULL;
  free(offsets); offsets = NULL;
  free(domainOffsets); domainOffsets = NULL;
  MPI_Finalize();

  return 0;
}
