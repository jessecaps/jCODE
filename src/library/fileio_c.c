/* -------------------------------------------------
   All binary files are in C format
   Fortran does not have a binary format by default
   
   The following functions handle
   - open
   - write
   - read
   - close
   ------------------------------------------------ */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

/* Constants */
/*const static int size_file = 64;*/
const static int none = 0;
const static int little = 1;
const static int big = 2;

/* Global variables */
static FILE* fp[100];
static int used[100] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

/* Folder list structure */
typedef struct element element;
struct element {
  char name[64];
  struct element * nxt;
};
typedef element* llist;
llist dirlist = NULL;

/* Definition of the functions */
void BINARY_FILE_OPEN(int* unit, const char* file, const char* mode, int* ierr, const int size1, const int size2);
void BINARY_FILE_CLOSE(int* unit, int* ierr);
void BINARY_FILE_REWIND(int* unit, int* ierr);
void BINARY_FILE_WRITE(int* unit, const void* buffer, int* count, int* size, int* ierr);
void BINARY_FILE_READ(int* unit, void* buffer, int* count, int* size, int* ierr);
void BINARY_BIGFILE_WRITE(int* unit, const void* buffer, long* count, int* size, int* ierr);
void BINARY_BIGFILE_READ(int* unit, void* buffer, long* count, int* size, int* ierr);
void CREATE_FOLDER(const char* file, const int size);
int FIND_IN_LIST(llist list, const char* name);
llist ADD_ELEMENT(llist list, const char* name);
void PRINT_LIST(llist list);

/* Fortran Wrapper*/
void binary_file_open(int* unit, const char* file, const char* mode, int* ierr, const int size1, const int size2){
  BINARY_FILE_OPEN(unit, file, mode, ierr, size1, size2);
}
void binary_file_open_(int* unit, const char* file, const char* mode, int* ierr, const int size1, const int size2){
  BINARY_FILE_OPEN(unit, file, mode, ierr, size1, size2);
}
void binary_file_open__(int* unit, const char* file, const char* mode, int* ierr, const int size1, const int size2){
  BINARY_FILE_OPEN(unit, file, mode, ierr, size1, size2);
}
void binary_file_close(int* unit, int* ierr){
  BINARY_FILE_CLOSE(unit, ierr);
}
void binary_file_close_(int* unit, int* ierr){
  BINARY_FILE_CLOSE(unit, ierr);
}
void binary_file_close__(int* unit, int* ierr){
  BINARY_FILE_CLOSE(unit, ierr);
}
void binary_file_rewind(int* unit, int* ierr){
  BINARY_FILE_REWIND(unit, ierr);
}
void binary_file_rewind_(int* unit, int* ierr){
  BINARY_FILE_REWIND(unit, ierr);
}
void binary_file_rewind__(int* unit, int* ierr){
  BINARY_FILE_REWIND(unit, ierr);
}
void binary_file_write(int* unit, const void* buffer, int* count, int* size, int* ierr){
  BINARY_FILE_WRITE(unit, buffer, count, size, ierr);
}
void binary_file_write_(int* unit, const void* buffer, int* count, int* size, int* ierr){
  BINARY_FILE_WRITE(unit, buffer, count, size, ierr);
}
void binary_file_write__(int* unit, const void* buffer, int* count, int* size, int* ierr){
  BINARY_FILE_WRITE(unit, buffer, count, size, ierr);
}
void binary_file_read(int* unit, void* buffer, int* count, int* size, int* ierr){
  BINARY_FILE_READ(unit, buffer, count, size, ierr);
}
void binary_file_read_(int* unit, void* buffer, int* count, int* size, int* ierr){
  BINARY_FILE_READ(unit, buffer, count, size, ierr);
}
void binary_file_read__(int* unit, void* buffer, int* count, int* size, int* ierr){
  BINARY_FILE_READ(unit, buffer, count, size, ierr);
}
void binary_bigfile_write(int* unit, const void* buffer, long* count, int* size, int* ierr){
  BINARY_BIGFILE_WRITE(unit, buffer, count, size, ierr);
}
void binary_bigfile_write_(int* unit, const void* buffer, long* count, int* size, int* ierr){
  BINARY_BIGFILE_WRITE(unit, buffer, count, size, ierr);
}
void binary_bigfile_write__(int* unit, const void* buffer, long* count, int* size, int* ierr){
  BINARY_BIGFILE_WRITE(unit, buffer, count, size, ierr);
}
void binary_bigfile_read(int* unit, void* buffer, long* count, int* size, int* ierr){
  BINARY_BIGFILE_READ(unit, buffer, count, size, ierr);
}
void binary_bigfile_read_(int* unit, void* buffer, long* count, int* size, int* ierr){
  BINARY_BIGFILE_READ(unit, buffer, count, size, ierr);
}
void binary_bigfile_read__(int* unit, void* buffer, long* count, int* size, int* ierr){
  BINARY_BIGFILE_READ(unit, buffer, count, size, ierr);
}
void create_folder(const char* name, const int size){
  CREATE_FOLDER(name,size);
}
void create_folder_(const char* name, const int size){
  CREATE_FOLDER(name,size);
}
void create_folder__(const char* name, const int size){
  CREATE_FOLDER(name,size);
}


/* Open the file */
void BINARY_FILE_OPEN(int* unit, const char* file, const char* mode, int* ierr, const int size1, const int size2){
  
  char *tmp;
  int test;
  
  /* Get file unit */
  for( *unit=0; *unit<100 && used[*unit]==1; (*unit)++) {}
  if (*unit==100) {
    printf("BINARY_FILE_OPEN: Trying to open too many files simultaneously\n");
    *ierr = 1;
    return;
  }
  used[*unit] = 1;

  /* Convert name from Fortran to C */
  char filename[64];
  strncpy(filename,file,(size_t) size1);
  filename[size1] = '\0';
  
  /* Check format and force binary */
  char format[3];
  if (mode[0]=='w')
    format[0] = 'w';
  else if (mode[0]=='r')
    format[0] = 'r';
  else if (mode[0]=='a')
    format[0] = 'a';
  else{
    printf("In fileio_c.c: BINARY_FILE_OPEN received an unknown format\n");
    exit(1);
  }
  format[1] = 'b';
  format[2] = '\0';
  
  /* Open file */
  fp[*unit] = fopen(filename,format);
  
  if (fp[*unit]==NULL){
    printf("BINARY_FILE_OPEN: Error opening the file\n");    
    *ierr = 2;
    return;
  }else
    *ierr = 0;
  
  return;
}

/* Close the file */
void BINARY_FILE_CLOSE(int* unit, int* ierr){
  
  fclose(fp[*unit]);
  used[*unit] = 0;
  return;
}

/* Rewind the file */
void BINARY_FILE_REWIND(int* unit, int* ierr){
  
  rewind(fp[*unit]);
  return;
}

/* Write to a file */
void BINARY_FILE_WRITE(int* unit, const void* buffer, int* count, int* size,  int* ierr){
  
  int i,n;
  unsigned char* tmp;
  char* pos;
  
  tmp = malloc(*size);
  pos = (char*) buffer;
  
  for (n=0; n<*count; n++){
    memcpy(tmp,pos,*size);
    for (i=0; i<*size; i++){
      fputc(tmp[i], fp[*unit]);
    }
    pos += *size;
  }
  
  free(tmp);
  *ierr = 0;
  return;
}

/* Read from a file */
void BINARY_FILE_READ(int* unit, void* buffer, int* count, int* size,  int* ierr){
  
  int i,n;
  unsigned char* tmp;
  char* pos;
  
  tmp = malloc(*size);
  pos = (char*) buffer;
  
  for (n=0; n<*count; n++){
    for (i=0; i<*size; i++){
      tmp[i] = fgetc(fp[*unit]);
    }
    memcpy(pos,tmp,*size);
    pos += *size;
  }
  
  free(tmp);
  *ierr = 0;
  return;
}

/* Write to a big file */
void BINARY_BIGFILE_WRITE(int* unit, const void* buffer, long* count, int* size,  int* ierr){
  
  int i;
  long n;
  unsigned char* tmp;
  char* pos;
  
  tmp = malloc(*size);
  pos = (char*) buffer;
  
  for (n=0; n<*count; n++){
    memcpy(tmp,pos,*size);
    for (i=0; i<*size; i++){
      fputc(tmp[i], fp[*unit]);
    }
    pos += *size;
  }
  
  free(tmp);
  *ierr = 0;
  return;
}

/* Read from a big file */
void BINARY_BIGFILE_READ(int* unit, void* buffer, long* count, int* size,  int* ierr){
  
  int i;
  long n;
  unsigned char* tmp;
  char* pos;
  
  tmp = malloc(*size);
  pos = (char*) buffer;
  
  for (n=0; n<*count; n++){
    for (i=0; i<*size; i++){
      tmp[i] = fgetc(fp[*unit]);
    }
    memcpy(pos,tmp,*size);
    pos += *size;
  }
  
  free(tmp);
  *ierr = 0;
  return;
}

/* Create a directory */
void CREATE_FOLDER(const char* name, const int size)
{
  int isThere,inList,isDir,status;
  struct stat buf;
  char dirname[64];
  
  /* process directory name */
  strncpy(dirname,name,(size_t) size);
  dirname[size]='\0';
  
  /* Check if directory "name" has been included in list */
  inList = FIND_IN_LIST(dirlist, dirname);
  if (inList == 1) {
    /*printf("%s: dir already found\n", dirname);*/
    return;
  }
  
  /* Check if directory "name" exists already */
  isThere = stat(dirname,&buf);
  isDir = S_ISDIR(buf.st_mode);
  if (isThere == 0 && isDir == 1){
    /*printf("%s: found directory, adding it to list\n",dirname);*/
    dirlist = ADD_ELEMENT(dirlist,dirname);
  }else if (isThere == 0 && isDir != 1){
    printf("%s: file and directory have same name!\n",dirname);
    exit(1);
  }else{
    dirlist = ADD_ELEMENT(dirlist,dirname);
#ifdef WIN32
    status = _mkdir(dirname);
#else
    status = mkdir(dirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
  }
  return;
}

/* Find an element in list */
int FIND_IN_LIST(llist list, const char* name)
{
  element *tmp = list;
  while(tmp != NULL) {
    if (strcmp(name,tmp->name) == 0){
      return 1;
    }else{
      tmp = tmp->nxt;
    }
  }
  return 0;
}

/* Add element to list */
llist ADD_ELEMENT(llist list, const char* name)
{  
  element* newElement = malloc(sizeof(element));
  strncpy(newElement->name,name,64);
  newElement->nxt = list;
  return newElement;
}

/* Print the list */
void PRINT_LIST(llist list)
{
  element *tmp = list;
  printf("List of stored dirnames:\n");
  while(tmp != NULL)
    {
      printf("check: %s\n", tmp->name);
      tmp = tmp->nxt;
    }
}
