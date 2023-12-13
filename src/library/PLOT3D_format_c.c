#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#define INT32_SWAPPED(x) ((((x)>>24)&0xff) |		\
			  (((x)>>8)&0xff00) |		\
			  (((x)<<8)&0xff0000) |		\
			  (((x)<<24)&0xff000000))
#define INT64_SWAPPED(x) ((((x)>>56)&0xff) |			\
			  (((x)>>40)&0xff00) |			\
			  (((x)>>24)&0xff0000) |		\
			  (((x)>>8)&0xff000000) |		\
			  (((x)<<8)&0xff00000000) |		\
			  (((x)<<24)&0xff0000000000) |		\
			  (((x)<<40)&0xff000000000000) |	\
			  (((x)<<56)&0xff00000000000000))

#define STRING_LENGTH 64
#define SIZEOF_SCALAR 8
#define SIZEOF_PLOT3D_OFF 8

#define plot3d_off_t int64_t
#define PLOT3D_OFF_SWAPPED(a) INT64_SWAPPED(a)

int detect_format(char filename[STRING_LENGTH + 1], int includeFunctionFiles, 
                  int32_t *nGrids, int *nDimensions, int *nScalars, int *isEndiannessNative, 
                  int *hasIblank, int *fileType, int **gridSizes);

int write_skeleton(char filename[STRING_LENGTH + 1], int32_t nGrids, int32_t nScalars, 
                   int fileType, int32_t *gridSizes);

void free_memory(int **gridSizes);

int detect_format_(FILE *file, int includeFunctionFiles, int32_t *nGrids, int *nScalars, 
                   int32_t **sizeHeader, int *isEndiannessNative, int *hasIblank, 
                   int *fileType);

int write_skeleton_(FILE *file, int32_t nGrids, int32_t nScalars, int fileType, 
                    int32_t *gridSizes);

int detect_format(char filename[STRING_LENGTH + 1], int includeFunctionFiles, int32_t *nGrids,
                  int *nDimensions, int *nScalars, int *isEndiannessNative, int *hasIblank, 
                  int *fileType, int **gridSizes) {

  int i, j, n, errorCode;
  int32_t *sizeHeader = NULL;
  FILE *file;
  int *gridSizes_;

  file = fopen(filename, "rb");
  if (!file) return -1; /* fopen failed */
  errorCode = detect_format_(file, includeFunctionFiles, nGrids, nScalars, 
                             &sizeHeader, isEndiannessNative, hasIblank, fileType);
  if (errorCode != 0) {
    if (sizeHeader != NULL)
      free(sizeHeader);    
    fclose(file);
    return errorCode;
  }

  nDimensions[0] = 1;
  n = (fileType[0] == 2) ? 4 : 3;
  for (i = 0; i < nGrids[0]; ++i)
    for (j = 0; j < 3; ++j)
      if (sizeHeader[j + n * i] > 1 && nDimensions[0] < j + 1) nDimensions[0] = j + 1;

  *gridSizes = (int*) malloc(nDimensions[0] * nGrids[0] * 4);
  gridSizes_ = *gridSizes;
  for (i = 0; i < nGrids[0]; ++i)
    for (j = 0; j < nDimensions[0]; ++j)
      gridSizes_[j + nDimensions[0] * i] = sizeHeader[j + n * i];

  free(sizeHeader);
  fclose(file);
  return 0;

}

int detect_format_(FILE *file, int includeFunctionFiles, int32_t *nGrids, int *nScalars, 
                   int32_t **sizeHeader, int *isEndiannessNative, int *hasIblank, 
                   int *fileType) {

  int i;
  int32_t *record;
  plot3d_off_t recordSize, recordSize_;

  if (!file) return -1;

  isEndiannessNative[0] = 1; /* assume data in file has same byte 
                                ordering as machine endianness */
  if (!fread(&recordSize, SIZEOF_PLOT3D_OFF, 1, file)) return -2; /* file is empty */
  if (recordSize != 4) { /* first record must be 4 bytes for a multi-block 
                            whole-format PLOT3D file */
    isEndiannessNative[0] = 0;
    recordSize = PLOT3D_OFF_SWAPPED(recordSize); /* try swapping */
    if (recordSize != 4) return -3; /* invalid file */
  }

  /* Check if the record markers are consistent, and the last record is followed by EOF */
  fseeko(file, 0, SEEK_SET);
  while (1) {
    if (!fread(&recordSize, SIZEOF_PLOT3D_OFF, 1, file)) { /* leading record size */
      if (feof(file)) break;
      return -2;
    }
    if (!isEndiannessNative[0]) 
      recordSize = PLOT3D_OFF_SWAPPED(recordSize); /* swap endianness if required */
    fseeko(file, recordSize, SEEK_CUR); /* skip the data */
    if (!fread(&recordSize_, SIZEOF_PLOT3D_OFF, 1, file))
      return -2; /* trailing record size */      
    if (!isEndiannessNative[0]) 
      recordSize_ = PLOT3D_OFF_SWAPPED(recordSize_); /* swap endianness if required */
    if (recordSize != recordSize_) return -4; /* leading/trailing record 
                                                 sizes don't match */
  }

  fseeko(file, SIZEOF_PLOT3D_OFF, SEEK_SET);
  if (!fread(nGrids, 4, 1, file)) return -2; /* first record contains number of grids */
  if (!isEndiannessNative[0]) 
    nGrids[0] = INT32_SWAPPED(nGrids[0]); /* swap endianness if required */
  if (*nGrids <= 0) return -3; /* must have positive number of grids */

  fseeko(file, SIZEOF_PLOT3D_OFF, SEEK_CUR);
  if (!fread(&recordSize, SIZEOF_PLOT3D_OFF, 1, file)) return -2; /* second record size */
  if (!isEndiannessNative[0]) 
    recordSize = PLOT3D_OFF_SWAPPED(recordSize); /* swap endianness if required */
  if (recordSize != 12 * nGrids[0] &&
      (!includeFunctionFiles || 
       recordSize != 16 * nGrids[0])) return -3; /* second record size is invalid */
  fileType[0] = (recordSize == 16 * nGrids[0]) ? 2 : 0; /* recognize that this is a function
                                                           file, or default to a grid file */

  *sizeHeader = (int32_t*) malloc(recordSize);
  record = *sizeHeader;

  if (!fread(record, 4, recordSize / 4, file)) return -2;
  if (!isEndiannessNative[0])
    for (i = 0; i < recordSize / 4; ++i)
      record[i] = INT32_SWAPPED(record[i]); /* swap endianness if required */
  for (i = 0; i < recordSize / 4; ++i)
    if (record[i] <= 0) return -3; /* grid size and number of components must be positive */
  if (fileType[0] == 2) {
    nScalars[0] = record[3];
    for (i = 0; i < nGrids[0]; ++i)
      if (record[i * 4 + 3] != nScalars[0]) return -3; /* number of components must 
                                                          be same for all grids */
  }

  fseeko(file, SIZEOF_PLOT3D_OFF, SEEK_CUR);
  if (fileType[0] == 2) {
    for (i = 0; i < nGrids[0]; ++i) {
      if (!fread(&recordSize, SIZEOF_PLOT3D_OFF, 1, file)) return -2;
      if (!isEndiannessNative[0])
        recordSize = PLOT3D_OFF_SWAPPED(recordSize); /* swap endianness if required */
      /* Check record size at the head of each grid's data */
      if (recordSize != SIZEOF_SCALAR * (plot3d_off_t)record[i * 4] * 
          (plot3d_off_t)record[i * 4 + 1] * (plot3d_off_t)record[i * 4 + 2] * 
          (plot3d_off_t)record[i * 4 + 3]) return -3;
      fseeko(file, SIZEOF_PLOT3D_OFF + recordSize, SEEK_CUR);
    }
  }
  else {
    for (i = 0; i < nGrids[0]; ++i) {
      if (!fread(&recordSize, SIZEOF_PLOT3D_OFF, 1, file)) return -2;
      if (!isEndiannessNative[0]) 
        recordSize = PLOT3D_OFF_SWAPPED(recordSize); /* swap endianness if required */
      if (recordSize == 4 * SIZEOF_SCALAR) { /* recognize that this is a solution file */
	fileType[0] = 1;
	fseeko(file, 4 * SIZEOF_SCALAR + SIZEOF_PLOT3D_OFF, SEEK_CUR);
	if (!fread(&recordSize, SIZEOF_PLOT3D_OFF, 1, file)) return -2;
	if (!isEndiannessNative[0]) 
          recordSize = PLOT3D_OFF_SWAPPED(recordSize); /* swap endianness if required */
      }
      /* Check record size at the head of each grid's data */
      if ((fileType[0] == 1 && recordSize != 5 * SIZEOF_SCALAR *
           (plot3d_off_t)record[i * 3] * (plot3d_off_t)record[i * 3 + 1] * 
           (plot3d_off_t)record[i * 3 + 2]) || 
	  (fileType[0] == 0 && (recordSize != 3 * SIZEOF_SCALAR * 
                                (plot3d_off_t)record[i * 3] * 
                                (plot3d_off_t)record[i * 3 + 1] * 
                                (plot3d_off_t)record[i * 3 + 2] && 
				recordSize != (3 * SIZEOF_SCALAR + 4) * 
                                (plot3d_off_t)record[i * 3] * 
				(plot3d_off_t)record[i * 3 + 1] * 
                                (plot3d_off_t)record[i * 3 + 2]))) return -3;
      hasIblank[0] = 0;
      if (fileType[0] == 0 && 
          recordSize == (3 * SIZEOF_SCALAR + 4) * 
          (plot3d_off_t)record[i * 3] * (plot3d_off_t)record[i * 3 + 1] * 
          (plot3d_off_t)record[i * 3 + 2]) hasIblank[0] = 1;
      fseeko(file, SIZEOF_PLOT3D_OFF + recordSize, SEEK_CUR);
    }
  }

  return 0;

}

int write_skeleton(char filename[STRING_LENGTH + 1], int32_t nGrids, int32_t nScalars, 
                   int fileType, int32_t *gridSizes) {

  int i, errorCode;
  FILE *file;
  
  file = fopen(filename, "wb");
  if (!file) return -1; /* fopen failed */
  errorCode = write_skeleton_(file, nGrids, nScalars, fileType, gridSizes);
  if (errorCode != 0) {
    fclose(file);
    return errorCode;
  }
  fclose(file);

  return 0;
  
}

int write_skeleton_(FILE *file, int32_t nGrids, int32_t nScalars, int fileType, 
                    int32_t *gridSizes) {

  int i, j;
  plot3d_off_t recordSize;

  if (!file) return -1;

  /* Write number of grids in Fortran unformatted style */
  recordSize = 4;
  if (fwrite(&recordSize, SIZEOF_PLOT3D_OFF, 1, file) != 1) return -2;
  if (fwrite(&nGrids, 4, 1, file) != 1) return -2;
  if (fwrite(&recordSize, SIZEOF_PLOT3D_OFF, 1, file) != 1) return -2;

  /* Grid size and number of components if this is a function file */
  recordSize = (fileType == 2) ? 16 * nGrids : 12 * nGrids;
  if (fwrite(&recordSize, SIZEOF_PLOT3D_OFF, 1, file) != 1) return -2;
  if (fileType == 2) {
    for (i = 0; i < nGrids; ++i) {
      for (j = 0; j < 3; ++j)
  	if (fwrite(&gridSizes[j + 3 * i], 4, 1, file) != 1) return -2;
      if (fwrite(&nScalars, 4, 1, file) != 1) return -2;
    }
  }
  else if (fwrite(gridSizes, 4, recordSize / 4, file) != recordSize / 4)
    return -2;
  if (fwrite(&recordSize, SIZEOF_PLOT3D_OFF, 1, file) != 1) return -2;

  if (fileType == 0) { /* grid file */
    for (i = 0; i < nGrids; ++i) {
      recordSize = (3 * SIZEOF_SCALAR + 4) * (plot3d_off_t)gridSizes[3 * i] * 
          (plot3d_off_t)gridSizes[3 * i + 1] * (plot3d_off_t)gridSizes[3 * i + 2];
      if (fwrite(&recordSize, SIZEOF_PLOT3D_OFF, 1, file) != 1) return -2;
      fseeko(file, recordSize, SEEK_CUR);
      if (fwrite(&recordSize, SIZEOF_PLOT3D_OFF, 1, file) != 1) return -2;
    }
  }
  else if (fileType == 1) { /* solution file */
    for (i = 0; i < nGrids; ++i) {
      recordSize = 4 * SIZEOF_SCALAR;
      if (fwrite(&recordSize, SIZEOF_PLOT3D_OFF, 1, file) != 1) return -2;
      fseeko(file, recordSize, SEEK_CUR);
      if (fwrite(&recordSize, SIZEOF_PLOT3D_OFF, 1, file) != 1) return -2;
      recordSize = 5 * SIZEOF_SCALAR * (plot3d_off_t)gridSizes[3 * i] * 
          (plot3d_off_t)gridSizes[3 * i + 1] * (plot3d_off_t)gridSizes[3 * i + 2];
      if (fwrite(&recordSize, SIZEOF_PLOT3D_OFF, 1, file) != 1) return -2;
      fseeko(file, recordSize, SEEK_CUR);
      if (fwrite(&recordSize, SIZEOF_PLOT3D_OFF, 1, file) != 1) return -2;
    }
  }
  else if (fileType == 2) { /* function file */
    for (i = 0; i < nGrids; ++i) {
      recordSize = (plot3d_off_t)nScalars * SIZEOF_SCALAR * (plot3d_off_t)gridSizes[3 * i] * 
          (plot3d_off_t)gridSizes[3 * i + 1] * (plot3d_off_t)gridSizes[3 * i + 2];
      if (fwrite(&recordSize, SIZEOF_PLOT3D_OFF, 1, file) != 1) return -2;
      fseeko(file, recordSize, SEEK_CUR);
      if (fwrite(&recordSize, SIZEOF_PLOT3D_OFF, 1, file) != 1) return -2;
    }
  }

  return 0;

}

void free_memory(int **gridSizes) {
  if (*gridSizes != NULL)
    free(*gridSizes);
}
