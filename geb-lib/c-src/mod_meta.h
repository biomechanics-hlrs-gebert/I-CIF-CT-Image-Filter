#ifndef MOD_META_H_
#define MOD_META_H_
/**
* \file
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
* \brief Header for the meta file format
*/

// ==== TODO's ====
/**
*   - maybe replace the macros with static gloabl variables
*   - implement use of different bitlengths for variables, maybe
*/

// ==== include statements ====
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>

// ==== macro definitions ====
// ---- file endings ----
#define LOG_SUFFIX ".log"
#define LOCK_SUFFIX ".lock"
#define HEAD_SUFFIX ".head"
#define META_SUFFIX ".meta"
#define MON_SUFFIX ".mon"
#define RES_SUFFIX ".result"
#define CSV_SUFFIX ".csv"
#define TEX_SUFFIX ".tex"
#define VTK_SUFFIX ".vtk"
#define RAW_SUFFIX ".raw"
// ---- META control macros ----
//datatypes.
#define META_MIK 4
#define META_IK 8
#define META_RK 8
//string lengths
#define META_MCL 512                    //meta max line length
#define META_SCL 64                     //meta short string length
#define META_KCL 25                     //Keyword max length
#define META_STDSPC 39                  //Keyword standard space (parameter max length)
#define META_UCL 8                      //Unit max length
//string symbols
#define META_LOCK_SUFFIX ".lock"
#define META_SUFFIX_DECLARATOR "."
#define META_BASENAME_SEPARATOR "_"
#define META_KEYWORD_SEPARATOR " "
#define META_SECTION_DECLARATOR "**"
#define META_KEYWORD_DECLARATOR "*"
//static array lengths
//TODO replace with dynamic allocation on the heap
#define META_MAX_FILE_LINES 1024
#define META_MAX_NUMBER_OF_KEYWORD_ITERATIONS 16
/**
* Define the behavior of multile occurences of meta file keyword
* Options are: 
*   - own (use the keyword in the programs own section)
*   - first(use the first occurrence)
*   - last (use the last occurrence)
*   - before (use keyword if it is defined in the sections before the own section, going back until its found)
*   - after (use keyword if it is defined in the sections after the own section, going forward until its found)
*
* syntax is: "PRIMARY[,SECONDARY]" where PRIMARY is the way to go and SECONDARY is the optional fallback option
*/
#define META_KEYWORD_OCCURRENCE_BEHAVIOUR "own,before"
/**
* If dimensions is not defined, define it.
*/
#ifndef DIMENSIONS
    #define DIMENSIONS 3
#endif
/**
* The size a meta file line may be long at maximum. It is composed of the individual specifiers
*/
#ifndef META_MCL
    #define META_MCL (META_KCL + META_STDSPC + META_UCL + META_TIME_MAX_LENGTH)
#else
    #if META_MCL < (META_KCL + META_STDSPC + META_UCL + META_TIME_MAX_LENGTH)
        #error META_MCL too small
    #endif
#endif

// ==== enum declarations ====
/**
* \typedef KeywordBehaviour;
* \enum KeywordBehaviour
*/
typedef enum {UNDEFINED, OWN, BEFORE, AFTER, FIRST, LAST} KeywordBehaviour;

// ==== struct declarations ====
/**
* \typedef basename
* \struct basename mod_meta.h "c-src/mod_meta.h"
*/
typedef struct {
    char full_name[META_MCL];
    char path[META_MCL];
    char path_and_basename[META_MCL];
    char basename[META_MCL];
    char dataset[META_MCL];
    char type[META_MCL];
    char purpose[META_MCL];
    char app[META_MCL];
    char features[META_MCL];
} basename;
/**
* \typedef metafile
* \struct metafile mod_meta.h "c-src/mod_meta.h"
*/
typedef struct {
    int number_of_lines;
    char content[META_MAX_FILE_LINES * META_MCL];
    basename *basename; //optional, add the basename to the metafile to bind these two together
} metafile;

// ==== global variables ====
/**
* \brief static global variables (only meta library internals)
*/
/* TODO
* use long variable names:
* WARNING: fh_meta_put was here renamed to fh_meta_out !!
*/
char * global_meta_program_keyword;
char * global_meta_prgrm_mstr_app;
static FILE * fh_meta_in;
static FILE * fh_meta_out;
static FILE * fh_mon;
static FILE * fh_out;
static FILE * fh_log;
static FILE * fh_res;
static FILE * fh_csv;
static FILE * fh_head;
static FILE * fh_tex;
static FILE * fh_vtk;
static FILE * fh_raw;
static basename in;
static basename out;
static clock_t meta_start;
static clock_t meta_end;


// ==== function declarations ====

//public interaction (read) function declarations
int meta_read_string(char *, metafile *, char *);
int meta_read_int_0D(char *, metafile *, int *);
int meta_read_int_1D(char *, metafile *, int dims, int[dims]);
int meta_read_double_0D(char *, metafile *, double *);
int meta_read_double_1D(char *, metafile *, int dims, double[dims]);

//public interactions (write) function declarations
int meta_write_int_0D(char *, char *, int);
int meta_write_int_1D(char *, char *, int dims, int[dims]);
int meta_write_long_0D(char *, char *, long long);
int meta_write_long_1D(char *, char *, int dims, long long[dims]);
int meta_write_double_0D(char *, char *, double);
int meta_write_double_1D(char *, char *, int dims, double[dims]);
int meta_write_string(char *, char *);

//public direct file interactions
int meta_handle_lock_file(char *, char *);
int meta_append(metafile *);
int meta_create_new(char *);
int meta_invoke(metafile *);
int meta_continue(metafile *);
int meta_start_ascii(FILE **, char *);
int meta_stop_ascii(FILE *, char *);
int meta_existing_ascii(FILE **, char *, int *);
int meta_signing(char *);
int meta_close(int);

//public helper function declarations
size_t meta_get_filesize(basename *);
size_t meta_count_lines(FILE *);
int meta_parse_basename(char *, char *);
int meta_check_unit(char *);
int meta_check_keyword(char *);
int meta_write_sha256sum(char *);
int meta_delete_empty_file(char *);
int meta_extract_keyword_data(char *, int dims, metafile *, char[dims][META_MCL]);
int meta_write_keyword(char *, char *, char *);

#endif
