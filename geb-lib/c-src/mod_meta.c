/**
* \file
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
* \brief Implementation file for the meta file format
*/

//include statements
#include "mod_meta.h"
#include "include_c/revision_meta.h"

//global variable definitions (declared in header)
/*global*/ basename in = {
    .full_name[0] = '\0',
    .path[0] = '\0',
    .path_and_basename[0] = '\0',
    .basename[0] = '\0',
    .dataset[0] = '\0',
    .type[0] = '\0',
    .purpose[0] = '\0',
    .app[0] = '\0',
    .features[0] = '\0'
};
/*global*/ basename out = {
    .full_name[0] = '\0',
    .path[0] = '\0',
    .path_and_basename[0] = '\0',
    .basename[0] = '\0',
    .dataset[0] = '\0',
    .type[0] = '\0',
    .purpose[0] = '\0',
    .app[0] = '\0',
    .features[0] = '\0'
};
/*global*/ char * global_meta_program_keyword = NULL;
/*global*/ char * global_meta_prgrm_mstr_app = NULL;

// ---- DECLARATIONS ----
static ssize_t __meta_get_filesize(char *);
static char *__meta_get_metafile_string_reference(metafile *, unsigned int);
static char *__meta_fast_strcat(char *, char *);
static void __meta_zero_basename_struct(basename *);
static void __meta_zero_array(size_t array_size, void *array);
static int __meta_get_file_suffix(char *filename_with_suffix, char *buffer);
static int __meta_fill_basename_struct(basename *);
static int __meta_replace_char(char *, char, char);
static int __meta_get_metafile_string_copy(metafile *, unsigned int, char *);
static int __meta_set_metafile_string(metafile *, unsigned int, char *);
static int __meta_print_error(FILE *, char *, int);

// ==== IMPLEMENTATIONS ====
// ---- public function implementations ----
/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method for parsing a character parameter from a keyword.
*
* \param[in] keyword The keyword to parse.
* \param[in] metafile A pointer to a metafile structure.
* \param[out] value Array where the value may be stored. Must be large enough.
* \return 0 on success, 1 otherwise.
*/
int meta_read_string(char *keyword, metafile *metafile, char *value){
    char tokens[1][META_MCL];
    
    if(meta_extract_keyword_data(keyword, 1, metafile, tokens)) 
        return 1;
    strcpy(value, tokens[0]);
    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method for parsing a integer parameter from a keyword.
*
* \param[in] keyword The keyword to parse.
* \param[in] metafile A pointer to a metafile structure.
* \param[out] value Pointer to integer where the value may be stored.
* \return 0 on success, 1 otherwise.
*/
int meta_read_int_0D(char *keyword, metafile *metafile, int *value){
    char tokens[1][META_MCL];
    
    if(meta_extract_keyword_data(keyword, 1, metafile, tokens))
        return 1;
    *value = atoi(tokens[0]);
    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method for parsing an array of integer parameters from a keyword.
*
* \param[in] keyword The keyword to parse.
* \param[in] metafile A pointer to a metafile structure.
* \param[in] dims Integer holding the number of array elements to be read.
* \param[out] values An array where the integers may be stored.
* \return 0 on success, 1 otherwise.
*/
int meta_read_int_1D(char *keyword, metafile *metafile, int dims, int values[dims]){
    char tokens[dims][META_MCL];

    if(meta_extract_keyword_data(keyword, dims, metafile, tokens))
        return 1;
    for(int i = 0; i < dims; i++)
        values[i] = atoi(tokens[i]);
    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method for parsing a floating point parameter from a keyword.
*
* \param[in] keyword The keyword to parse.
* \param[in] metafile A pointer to a metafile structure.
* \param[out] value Pointer to a float where the value may be stored.
* \return 0 on success, 1 otherwise.
*/
int meta_read_double_0D(char *keyword, metafile *metafile, double *value){
    char tokens[1][META_MCL];
    
    if(meta_extract_keyword_data(keyword, 1, metafile, tokens))
        return 1;
    *value = atof(tokens[0]);

    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method for parsing an array on array of floating point parameters from a keyword.
*
* \param[in] keyword The keyword to parse.
* \param[in] metafile A pointer to a metafile structure.
* \param[in] dims Integer holding the number of array elements to be read.
* \param[out] values Array where the values may be stored.
* \return 0 on success, 1 otherwise.
*/
int meta_read_double_1D(char *keyword, metafile *metafile, int dims, double values[dims]){
    char tokens[dims][META_MCL];
    
    if(meta_extract_keyword_data(keyword, dims, metafile, tokens)) 
        return 1;
    for(int i = 0; i < dims; i++)
        values[i] = atof(tokens[i]);
    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method for writing a keyword with corresponding integer scalar value.
*
* \param[in] keyword Keyword to write to the metafile.
* \param[in] unit Unit to write behind the data. Nullable.
* \param[in] value The integer scalar value to be written behind the keyword.
* \return 0 on success, 1 otherwise.
*/
int meta_write_int_0D(char *keyword, char *unit, int value){
    char buffer[META_STDSPC];

    if(snprintf(buffer, META_STDSPC, "%d", value) == META_STDSPC) 
        return 1;
    return meta_write_keyword(keyword, buffer, unit == NULL ? "" : unit);
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method for writing a keyword with corresponding integer vector value.
*
* \param[in] keyword Keyword to write to the metafile.
* \param[in] unit Unit to write behind the data. Nullable.
* \param[in] dims Integer holding the number of arguments to be written.
* \param[in] value The integer vector value to be written behind the keyword.
* \return 0 on success, 1 otherwise.
*/
int meta_write_int_1D(char *keyword, char *unit, int dims, int value[dims]){
    if(value == NULL)
        return 1;

    char buffer[META_STDSPC], tmp[META_STDSPC];
    size_t string_length = 0;
    buffer[0] = '\0';
    for(int i = 0; i < dims; i++){
        if((string_length += snprintf(tmp, META_STDSPC, "%d ", value[i])) >= META_STDSPC)
            return 1;
        strcat(buffer, tmp);
    }
    return meta_write_keyword(keyword, buffer, unit == NULL ? "" : unit);
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method for writing a keyword with corresponding long scalar value.
*
* \param[in] keyword Keyword to write to the metafile.
* \param[in] unit Unit to write behind the data. Nullable.
* \param[in] value The long scalar value to be written behind the keyword.
* \return 0 on success, 1 otherwise.
*/
int meta_write_long_0D(char *keyword, char *unit, long long value){
    char buffer[META_STDSPC];
    if(snprintf(buffer, META_STDSPC, "%lld", value) == META_STDSPC)
        return 1;
    return meta_write_keyword(keyword, buffer, unit == NULL ? "" : unit);
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method for writing a keyword with corresponding long vector value.
*
* \param[in] keyword Keyword to write to the metafile.
* \param[in] unit Unit to write behind the data. Nullable.
* \param[in] dims Integer holding the number of arguments to be written.
* \param[in] value The long vector value to be written behind the keyword.
* \return 0 on success, 1 otherwise.
*/
int meta_write_long_1D(char *keyword, char *unit, int dims, long long value[dims]){
    if(value == NULL) 
        return 1;

    char buffer[META_STDSPC], tmp[META_STDSPC];
    size_t string_length = 0;
    buffer[0] = '\0';
    for(int i = 0; i < dims; i++){
        if((string_length += snprintf(tmp, META_STDSPC, "%lld ", value[i])) >= META_STDSPC)
            return 1;
        strcat(buffer, tmp);
    }
    return meta_write_keyword(keyword, buffer, unit == NULL ? "" : unit);
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method for writing a keyword with corresponding double scalar value.
*
* \param[in] keyword Keyword to write to the metafile.
* \param[in] unit Unit to write behind the data. Nullable.
* \param[in] value The double scalar value to be written behind the keyword.
* \return 0 on success, 1 otherwise.
*/
int meta_write_double_0D(char *keyword, char *unit, double value){
    char buffer[META_STDSPC];
    if(snprintf(buffer, META_STDSPC, "%lf", value) == META_STDSPC)
        return 1;
    return meta_write_keyword(keyword, buffer, unit == NULL ? "" : unit);
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method for writing a keyword with corresponding double vector value.
*
* \param[in] keyword Keyword to write to the metafile.
* \param[in] unit Unit to write behind the data. Nullable.
* \param[in] dims Integer holding the number of arguments to be written.
* \param[in] value The double vector value to be written behind the keyword.
* \return 0 on success, 1 otherwise.
*/
int meta_write_double_1D(char *keyword, char *unit, int dims, double value[dims]){
    if(value == NULL)
        return 1;
    char buffer[META_STDSPC], tmp[META_STDSPC];
    size_t string_length = 0;
    buffer[0] = '\0';
    for(int i = 0; i < dims; i++){
        if((string_length += snprintf(tmp, META_STDSPC, "%lf ", value[i])) >= META_STDSPC)
            return 1;
        strcat(buffer, tmp);
    }
    return meta_write_keyword(keyword, buffer, unit == NULL ? "" : unit);
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method for writing a keyword with corresponding string value.
*
* \param[in] keyword Keyword to write to the metafile.
* \param[in] unit Unit to write behind the data. Nullable.
* \param[in] value The string value to be written behind the keyword.
* \return 0 on success, 1 otherwise.
*/
int meta_write_string(char *keyword, char *value){
    if(value == NULL)
        return 1;
    if(strnlen(value, META_STDSPC) == META_STDSPC)
        return 1;
    return meta_write_keyword(keyword, value, "");
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method to encapsulate the lock file handling.
*
* Descritpion:
*   The lock file might be used as an activity tracker for large
*   computations.
*
* \param[in] restart Wether to restart or not to.
* \param[in] restart_cmdarg Nullable. Possible cmd argument override.
* \return 0 on success, 1 otherwise.
*/
int meta_handle_lock_file(char *restart, char *restart_cmdarg){
    if(restart == NULL)
        return 1;
    if(LOCK_SUFFIX == NULL)
        return 1;
    
    bool exist = false;
    int error;
    char lockname[2 * META_MCL + strlen(LOCK_SUFFIX) + 1], command[META_MCL + 3], *ptr;
    char *message;
    
    if(restart_cmdarg != NULL){
        if(meta_compare_restart(restart, restart_cmdarg))
            return 1;
    }

    ptr = lockname;
    ptr = __meta_fast_strcat(ptr, /*global*/ in.path);
    ptr = __meta_fast_strcat(ptr, ".");
    ptr = __meta_fast_strcat(ptr, /*global*/ in.basename);
    ptr = __meta_fast_strcat(ptr, LOCK_SUFFIX);

    if(ptr == NULL)
        return 1;
    if(__meta_get_filesize(lockname) != -1 && !strcmp(restart, "N")){
        if(__meta_get_filesize(/*global*/ out.full_name) != -1){
            ptr = __meta_fast_strcat(command, "rm ");
            __meta_fast_strcat(ptr, /*global*/ out.full_name);
            system(command);
        }
        return __meta_print_error(stdout, "The .*.lock file is set and a restart is prohibited by default or the user.\n", 1);
    }

    if(__meta_get_filesize(lockname) == -1 && (!strcmp(restart, "Y") || !strcmp(restart, "N"))){
        ptr = __meta_fast_strcat(command, "touch ");
        __meta_fast_strcat(ptr, lockname);
        if(error = system(command))
            return __meta_print_error(stdout, "The .*.lock file could not be set.\n", error);
    }
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Compare restart arguments.
*
* \param[inout] restart Wether to restart or not to.
* \param[in] restart_cmdarg Nullable. Possible cmd argument override.
* \return 0 on success, 1 otherwise.
*/
int meta_compare_restart(char *restart, char *restart_cmdarg){
        if(restart == NULL || restart_cmdarg == NULL)
            return 1;
        
        char *message;

        if(strcmp(restart_cmdarg, "") && strcmp(restart_cmdarg, "U")){
            if(!strcmp(restart_cmdarg, "N")){
                strcpy(restart, restart_cmdarg);
                message = "The keyword »restart« was overwritten by the command flag --no-restart\n";
            }
            else if(!strcmp(restart_cmdarg, "Y")){
                message = "The keyword »restart« was overwritten by the command flag --restart\n";
                strcpy(restart, restart_cmdarg);
            }
            else
                message = "";
                fprintf(stdout, message);
        }
        return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Module to extract the data string of keywords.
*
* Description:
*   Module to parse information of keywords.
*   An arbitrary keyword with up to META_KCL characters may be specified.
*   The program reads keywords as long as they are before or within the owns
*   program scope.
*   
* \param[in] keyword Keyword to read.
* \param[in] dims Dimension requested.
* \param[in] m_in Metafile structure to read keyword from.
* \param[out] res_tokens Results stored in array of dimension META_MCL times dims.
* \return 0 on success, 1 otherwise
*/
int meta_extract_keyword_data(char *keyword, int dims, metafile *m_in, char res_tokens[dims][META_MCL]){
    if(keyword == NULL || m_in == NULL || res_tokens == NULL) 
        return 1;
    if(/*global*/ global_meta_program_keyword == NULL)
        return 1;
    if(dims < 1)
        return 1;

    size_t i, j;
    bool own_section = false, keyword_found = false;
    char buffer[META_MCL], kwargs_storage[dims][META_STDSPC], *token;
    
    meta_check_keyword(keyword);
    __meta_zero_array(dims * META_STDSPC, (void *) kwargs_storage);

    for(i = 0; i < m_in -> number_of_lines; i++){
        if(__meta_get_metafile_string_copy(m_in, i, buffer))
            return 1; 
        token = strtok(buffer, META_KEYWORD_SEPARATOR);
        if(token == NULL)
            return 1;
        if(token[0] == '\n' || token[0] == '\0') 
            continue; 

        //handle new section
        if(!strcmp(token, META_SECTION_DECLARATOR)){
            token = strtok(NULL, META_KEYWORD_SEPARATOR);
            __meta_replace_char(token, '\n', '\0');
            if(own_section)
                break;
            if(!strcmp(token, /*global*/ global_meta_program_keyword)) 
                own_section = true;
            continue;
        }

        //handle new keyword
        else if(!strcmp(token, META_KEYWORD_READ_DECLARATOR) || !strcmp(token, META_KEYWORD_WRITE_DECLARATOR)){
            token = strtok(NULL, META_KEYWORD_SEPARATOR);
            if(token == NULL) 
                return 1; //Error: no keyword specified
            if(!strcmp(token, keyword)){
                token = strtok(NULL, META_KEYWORD_SEPARATOR);
                keyword_found = true;
                for(j = 0; j < dims; j++){
                    if(token == NULL) 
                        return 1; //error: there are not enough values to describe all dims
                    if(strnlen(token, META_STDSPC) == META_STDSPC)
                        return 1; //keyword parameter too long
                    strcpy(kwargs_storage[j], token);
                    token = strtok(NULL, META_KEYWORD_SEPARATOR);
                }
            }
        }
        else
            continue; //skip the uninterpretable line
    }
    
    if(!keyword_found) 
        return 1; //the whole file does not contain the required keyword.
    for(i = 0; i < dims; i++){
        __meta_replace_char(kwargs_storage[i], '\n', '\0');
        strcpy(res_tokens[i], kwargs_storage[i]);
    }
    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method to write finalized strings of keywords.
*
* \param[in] keyword Keyword to write.
* \param[in] stdspcfill String with data.
* \param[in] unit Unit of the value.
* \return 0 on success, 1 otherwise.
*/
int meta_write_keyword(char *keyword, char *stdspcfill, char *unit){
    if(keyword == NULL || stdspcfill == NULL || unit == NULL) 
        return 1;
    if(meta_check_keyword(keyword) < 0)
        return 1;
    if(meta_check_unit(unit) < 0)
        return 1;

    if(strnlen(stdspcfill, META_STDSPC) == META_STDSPC) 
        return 1;

    char metafile_line[META_MCL], time_buffer[META_MCL], *ptr, *max;
    time_t now_time;
    struct tm *time_struct;

    now_time = time(NULL);
    time_struct = localtime(&now_time);
    if(strftime(time_buffer, sizeof(time_buffer), "%c", time_struct) >= META_SCL) 
        return 1;
    
    ptr = metafile_line;
    ptr = __meta_fast_strcat(ptr, META_KEYWORD_WRITE_DECLARATOR);
    ptr = __meta_fast_strcat(ptr, META_KEYWORD_SEPARATOR);
    ptr = __meta_fast_strcat(ptr, keyword);
    for(max = ptr + META_KCL - strlen(keyword); ptr < max; ptr += strlen(META_KEYWORD_SEPARATOR) - 1)
        ptr = __meta_fast_strcat(ptr, META_KEYWORD_SEPARATOR);
    ptr = __meta_fast_strcat(ptr, META_KEYWORD_SEPARATOR);
    ptr = __meta_fast_strcat(ptr, stdspcfill);
    for(max = ptr + META_STDSPC - strlen(stdspcfill); ptr < max; ptr += strlen(META_KEYWORD_SEPARATOR) - 1)
        ptr = __meta_fast_strcat(ptr, META_KEYWORD_SEPARATOR);
    ptr = __meta_fast_strcat(ptr, META_KEYWORD_SEPARATOR);
    ptr = __meta_fast_strcat(ptr, unit);
    for(max = ptr + META_UCL - strlen(unit); ptr < max; ptr += strlen(META_KEYWORD_SEPARATOR) -1)
        ptr = __meta_fast_strcat(ptr, META_KEYWORD_SEPARATOR);
    ptr = __meta_fast_strcat(ptr, META_KEYWORD_SEPARATOR);
    ptr = __meta_fast_strcat(ptr, time_buffer);
    ptr = __meta_fast_strcat(ptr, "\n");

    if(fseek(/*global*/ fh_meta_out, 0, SEEK_END))
        return 1;
    if(((ptr - metafile_line) / sizeof(char)) != fprintf(/*global*/ fh_meta_out, metafile_line)) 
        return 1;
    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
* 
* \brief Method to open a meta file to append data / keywords.
* 
* \param[out] metafile A pointer to a metafile structure.
* \param[in] size_mpi The number of processors invloved in the computation.
* \return 0 on success, 1 on error.
*/
/*
* UNIFY HEAD:             |  c error(type=metafile)
* UNIFY DOC:              |  fortran error(inout=meta_as_rry), c error(type=metafile)
* FUNCTIONAL (THEORY):    |  done
* FUNCTIONAL (PRACTICE):  |  done
* UNIT TESTING:           |  
* ADD  ERROR HANDLER:     |  
* ADD  INPUT SANITIZER:   |  done
* MAKE PRETTY:            |  done
*/
int meta_append(metafile *metafile, int size_mpi){
    if(metafile == NULL)
        return 1;

    if(meta_invoke(metafile)) 
        return 1;
    if(meta_continue(metafile, size_mpi)) 
        return 1;
    if(meta_write_string("META_PARSED/INVOKED", "Now - Date/Time on the right."))
        return 1;

    /*global*/ meta_start = clock();

    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method to create a new meta file.
*
* \param[in] filename_with_suffix The filename for the new metafile.
* \return 0 on success, 1 otherwise.
*/
/*
* UNIFY HEAD:             |  done
* UNIFY DOC:              |  fortran error
* FUNCTIONAL (THEORY):    |  done
* FUNCTIONAL (PRACTICE):  |  
* UNIT TESTING:           |  
* ADD  ERROR HANDLER:     |  
* ADD  INPUT SANITIZER:   |  done
* MAKE PRETTY:            |  done
*/
int meta_create_new(char *filename_with_suffix){
    if(filename_with_suffix == NULL)
        return 1;

    size_t filename_length = strlen(filename_with_suffix);
    size_t meta_suffix_length = strlen(META_SUFFIX);
    size_t meta_out_length = strlen(/*global*/ out.full_name);
    char buffer[filename_length + meta_suffix_length + 1], *ptr;
    char extension_buffer[filename_length];
    char errormessage[
        filename_length >= meta_out_length + meta_suffix_length ?
        filename_length + 27 :
        meta_out_length + meta_suffix_length + 27
    ];

    if(__meta_get_filesize(filename_with_suffix) == -1){
        ptr = __meta_fast_strcat(errormessage, "The file ");
        ptr = __meta_fast_strcat(ptr, filename_with_suffix);
        __meta_fast_strcat(ptr, " does not exist.\n");
        return __meta_print_error(stdout, errormessage, 1);
    }
    
    if(__meta_get_file_suffix(filename_with_suffix, extension_buffer))
        return 1;
    if(meta_parse_basename(filename_with_suffix, extension_buffer)) 
        return 1;

    ptr = __meta_fast_strcat(buffer, /*global*/ in.path_and_basename);
    __meta_fast_strcat(ptr, META_SUFFIX);
    if(__meta_get_filesize(buffer) != -1){
        ptr = __meta_fast_strcat(errormessage, "The file ");
        ptr = __meta_fast_strcat(ptr, buffer);
        __meta_fast_strcat(ptr, " already exists.\n");
        return __meta_print_error(stdout, errormessage, 1);
    }

    /*global*/ fh_meta_out = fopen(buffer, "w+");
    if(/*global*/ fh_meta_out == NULL) 
        return 1;
    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method to open and prepare a meta file for use.
*
* \param[inout] metafile Pointer to a metafile structure.
* \return 0 on success, 1 otherwise.
*/
int meta_invoke(metafile *metafile){
    long long lines;
    if(metafile == NULL) 
        return 1;
    if(__meta_get_filesize(/*global*/ in.full_name) == -1)
        return 1;

    size_t buffer_size = 0;
    char meta_suffix[META_MCL];
    char filename_buffer[strlen(meta_suffix) + strlen(META_SUFFIX_DECLARATOR) + 1], *line_buffer = NULL;
    
    if(__meta_get_file_suffix(/*global*/ in.full_name, meta_suffix))
        return 1;
    strcpy(filename_buffer, META_SUFFIX_DECLARATOR);
    strcat(filename_buffer, meta_suffix);
    if(strcmp(filename_buffer, META_SUFFIX) != 0)
        return 1;
    if(meta_parse_basename(/*global*/ in.full_name, META_SUFFIX))
        return 1;

    /*global*/ fh_meta_in = fopen(/*global*/ in.full_name, "r+");
    
    lines = meta_count_lines(/*global*/ fh_meta_in);
    if(lines == -1)
        return 1;
    
    //allocate large enough buffer to hold all lines
    metafile -> content = (char *) calloc(lines * META_MCL, sizeof(char));
    metafile -> number_of_lines = lines;

    for(int i = 0; i < lines; i++){
        if(getline(&line_buffer, &buffer_size, /*global*/ fh_meta_in) == -1){
            free(line_buffer);
            return 1; 
        }
        if(buffer_size >= META_MCL){
            free(line_buffer);
            return 1;
        }
        if(__meta_set_metafile_string(metafile, i, line_buffer)){
            free(line_buffer);
            return 1;
        }
    }
    free(line_buffer);
    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method to invoke the meta output file.
*
* \param[inout] metafile Pointer to a metafile structure.
* \param[in] size_mpi The number of processors involved in the computation.
* \return 0 on success, 1 otherwise.
*/
int meta_continue(metafile *metafile, int size_mpi){
    if(metafile == NULL) 
        return 1;
    if(size_mpi < 1)
        return 1;

    int error;
    size_t stringlength;
    char *ptr;

    if(meta_read_string("NEW_BSNM_FEATURE", metafile, /*global*/ out.features)) 
        return 1;
    if(meta_read_string("NEW_BSNM_PURPOSE", metafile, /*global*/ out.purpose)) 
        return 1;

    if(
        !strcmp(/*global*/ out.purpose, /*global*/ in.purpose) && 
        !strcmp(/*global*/ out.features, /*global*/ in.features)
    )
        fprintf(stdout, "The basename (in part) did not change.");
    
    stringlength = strlen(/*global*/ out.path);
    stringlength += strlen(/*global*/ out.dataset);
    stringlength += strlen(/*global*/ out.type);
    stringlength += strlen(/*global*/ out.purpose);
    stringlength += strlen(/*global*/ out.features);
    stringlength += strlen(/*global*/ global_meta_prgrm_mstr_app);
    stringlength += 4 * strlen(META_BASENAME_SEPARATOR);
    if(stringlength >= META_MCL)
        return 1;

    ptr = __meta_fast_strcat(/*global*/ out.basename, /*global*/ out.dataset);
    ptr = __meta_fast_strcat(ptr, META_BASENAME_SEPARATOR);
    ptr = __meta_fast_strcat(ptr, /*global*/ out.type);
    ptr = __meta_fast_strcat(ptr, META_BASENAME_SEPARATOR);
    ptr = __meta_fast_strcat(ptr, /*global*/ out.purpose);
    ptr = __meta_fast_strcat(ptr, META_BASENAME_SEPARATOR);
    ptr = __meta_fast_strcat(ptr, /*global*/ global_meta_prgrm_mstr_app);
    ptr = __meta_fast_strcat(ptr, META_BASENAME_SEPARATOR);
    ptr = __meta_fast_strcat(ptr, /*global*/ out.features);
    if(ptr == NULL)
        return 1;

    ptr = __meta_fast_strcat(/*global*/ out.path_and_basename, /*global*/ out.path);
    ptr = __meta_fast_strcat(ptr, /*global*/ out.basename);
    if(ptr == NULL)
        return 1;
    
    ptr = __meta_fast_strcat(/*global*/ out.full_name, /*global*/ out.path_and_basename);
    ptr = __meta_fast_strcat(ptr, META_SUFFIX);
    if(ptr == NULL)
        return 1;

    char command[strlen(/*global*/ in.full_name) + strlen(/*global*/ out.full_name) + 5];
    strcpy(command, "cp ");
    strcat(command, /*global*/ in.full_name);
    strcat(command, " ");
    strcat(command, /*global*/ out.full_name);
    if(error = system(command))
        return __meta_print_error(stdout, "The update of the meta filename went wrong.\n", error);
    
    if(/*global*/ fh_meta_out != NULL) 
        return 1;
    /*global*/ fh_meta_out = fopen(out.full_name, "a");
    if(/*global*/ fh_meta_out == NULL) 
        return 1;

    if(fseek(/*global*/ fh_meta_out, 0, SEEK_END))
        return 1;

    if(meta_write_int_0D("PROCESSORS", NULL, size_mpi))
        return 1;

    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method to deal with the logging and renaming of the log in context 
*        of the meta data approach.
*
* Description:
*   This subroutine only gets called after meta append or meta close respectively.
*   If the variable in/out are not set, the program will not start / stop accordingly.
*   During computation (before meta_stop_ascii), the files are called 'temporary.suffix'
*   to show which ones are subject to modifications.
*    
* \param[out] fh File handle pointer to be filled with temporary file.
* \param[in] suf Suffix of the file.
* \return 0 on success, 1 otherwise.
*/
/*
* UNIFY HEAD:             |  done
* UNIFY DOC:              |  fortran error(inout_error=fh)
* FUNCTIONAL (THEORY):    |  done, c error(differ=true["error printing is false in fortran"])
* FUNCTIONAL (PRACTICE):  |  
* UNIT TESTING:           |  
* ADD  ERROR HANDLER:     |  done
* ADD  INPUT SANITIZER:   |  done
* MAKE PRETTY:            |  done
*/
int meta_start_ascii(FILE **fh, char *suf){
    if(fh == NULL || suf == NULL) 
        return 1;
    
    size_t temp_size = strlen(/*global*/ out.path) + strlen(suf) + 10;
    size_t perm_size = strlen(/*global*/ out.path_and_basename) + strlen(suf) + 1;
    char temporary_filename[temp_size], permanent_filename[perm_size];
    char temp_command[temp_size + 6], perm_command[perm_size + 6];
    char errortext[temp_size > perm_size ? temp_size + 21 : perm_size + 21];
    char *ptr;
    int error;

    strcpy(temporary_filename, /*global*/ out.path);
    strcat(temporary_filename, "temporary");
    strcat(temporary_filename, suf);

    strcpy(permanent_filename, /*global*/ out.path_and_basename);
    strcat(permanent_filename, suf);

    if(__meta_get_filesize(temporary_filename) != -1){
        strcpy(temp_command, "rm -r ");
        strcat(temp_command, temporary_filename);
        if(error = system(temp_command)){
            ptr = __meta_fast_strcat(errortext, "»");
            ptr = __meta_fast_strcat(ptr, temporary_filename);
            ptr = __meta_fast_strcat(ptr, "« not deletable.\n");
            return __meta_print_error(stdout, errortext, error);
        }
    }

    if(__meta_get_filesize(permanent_filename) != -1){
        strcpy(perm_command, "rm -r ");
        strcat(perm_command, permanent_filename);
        if(error = system(perm_command)){
            ptr = __meta_fast_strcat(errortext, "»");
            ptr = __meta_fast_strcat(ptr, permanent_filename);
            ptr = __meta_fast_strcat(ptr, "« not deletable.\n");
            return __meta_print_error(stdout, errortext, error);
        }
    }

    *fh = fopen(temporary_filename, "w");
    if(*fh == NULL)
        return 1;
    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method to stop the logging and renaming of additional ascii files.
*
* \param[in] fh File handle to the input file.
* \param[in] suf Suffix of the file.
* \return 0 on success, 1 otherwise.
*/
int meta_stop_ascii(FILE *fh, char *suf){
    if(fh == NULL || suf == NULL) 
        return 1;
    
    size_t temp_size = strlen(/*global*/ out.path) + strlen(suf) + 10;
    size_t perm_size = strlen(/*global*/ out.path_and_basename) + strlen(suf) + 1;
    size_t alt_size = strlen(/*global*/ in.path_and_basename) + strlen(suf) + 1;
    char temporary_filename[temp_size], permanent_filename[perm_size], alternative_filename[alt_size];
    char exist_command[temp_size + perm_size + 5], not_exist_command[alt_size + perm_size + 5];
    char errortext[temp_size + 68], *ptr;
    int error;

    ptr = __meta_fast_strcat(temporary_filename, /*global*/ out.path);
    ptr = __meta_fast_strcat(ptr, "temporary");
    __meta_fast_strcat(ptr, suf);

    ptr = __meta_fast_strcat(permanent_filename, /*global*/ out.path_and_basename);
    __meta_fast_strcat(ptr, suf);

    ptr = __meta_fast_strcat(alternative_filename, /*global*/ in.path_and_basename);
    __meta_fast_strcat(ptr, suf);

    if(fclose(fh) == EOF) 
        return 1;

    if(__meta_get_filesize(temporary_filename) != -1){
        ptr = __meta_fast_strcat(exist_command, "mv ");
        ptr = __meta_fast_strcat(ptr, temporary_filename);
        ptr = __meta_fast_strcat(ptr, " ");
        __meta_fast_strcat(ptr, permanent_filename);
        if(error = system(exist_command)){
            ptr = __meta_fast_strcat(errortext, "Can not rename the suffix_file from »");
            ptr = __meta_fast_strcat(ptr, temporary_filename);
            __meta_fast_strcat(ptr, "« to the proper basename.\n");
            return __meta_print_error(stdout, errortext, error);
        }
    }
    else{
        if(__meta_get_filesize(alternative_filename) != -1){
            ptr = __meta_fast_strcat(not_exist_command, "cp ");
            ptr = __meta_fast_strcat(ptr, alternative_filename);
            ptr = __meta_fast_strcat(ptr, " ");
            __meta_fast_strcat(ptr, permanent_filename);
            if(error = system(not_exist_command)){
                ptr = __meta_fast_strcat(errortext, "Can not copy the suffix_file from »");
                ptr = __meta_fast_strcat(ptr, temporary_filename);
                __meta_fast_strcat(ptr, "« to the proper basename.\n");
                return __meta_print_error(stdout, errortext, error);
            }
        }
    }
    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method to check and open ascii files wich must exist.
*
* Description:
*   For example to read input data. To stop the file, use meta_stop_ascii.
*
* \param[inout] fh File handle to the input file.
* \param[in] suf Suffix of the file.
* \param[out] amount_of_lines Integer pointer to store the number of lines in file.
* \return 0 on success, 1 otherwise.
*/
int meta_existing_ascii(FILE **fh, char *suf, int *amount_of_lines){
    if(fh == NULL || suf == NULL || amount_of_lines == NULL) 
        return 1;

    size_t filename_size = strlen(/*global*/ out.path) + strlen(suf) + 10;
    char filename[filename_size], errortext[filename_size + 43], *ptr;
    int amount_of_lines_temp;
    
    ptr = __meta_fast_strcat(filename, /*global*/ out.path);
    ptr = __meta_fast_strcat(ptr, "temporary");
    ptr = __meta_fast_strcat(ptr, suf);

    if(__meta_get_filesize(filename) == -1){
        ptr = __meta_fast_strcat(errortext, "The input file requested does not exist: ");
        ptr = __meta_fast_strcat(ptr, filename);
        ptr = __meta_fast_strcat(ptr, "\n");
        return __meta_print_error(stdout, errortext, 1);
    }

    *fh = fopen(filename, "r+");
    if(*fh == NULL) 
        return 1;
    amount_of_lines_temp = meta_count_lines(*fh);
    if(amount_of_lines_temp == -1) 
        return 1;
    *amount_of_lines = amount_of_lines_temp;
    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method to sign a meta file at the end of computation.
*
* \param[in] binary The file of the to be hashed binary.
* \return 0 on success, 1 otherwise.
*/
/*
* UNIFY HEAD:             |  done
* UNIFY DOC:              |  fortran error(unclear)
* FUNCTIONAL (THEORY):    |  done
* FUNCTIONAL (PRACTICE):  | 
* UNIT TESTING:           | 
* ADD  ERROR HANDLER:     |  unnessesary
* ADD  INPUT SANITIZER:   |  done
* MAKE PRETTY:            |  done
*/
int meta_signing(char *binary){
    if(binary == NULL) 
        return 1;
    
    size_t binary_length = strlen(binary);
    char *token, binary_copy[binary_length + 1];
    strcpy(binary_copy, binary);

    token = strtok(binary_copy, "_");
    if(token == NULL)
        return 1;
    token = strtok(NULL, "_");
    if(token == NULL)
        return 1;

    if(meta_write_string("PROGRAM_VERSION", token))
        return 1;
    if(meta_write_string("PROGRAM_GIT_HASH", /*global*/ hash))
        return 1;
    if(meta_write_sha256sum(binary))
        return 1;

    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method to close a meta file.
*
* \return 0 on success, 1 otherwise.
*/
/*
* UNIFY HEAD:             |  done
* UNIFY DOC:              |  c error(missing=description)
* FUNCTIONAL (THEORY):    |  done
* FUNCTIONAL (PRACTICE):  |  done
* UNIT TESTING:           | 
* ADD  ERROR HANDLER:     |  unnessesary
* ADD  INPUT SANITIZER:   |  done
* MAKE PRETTY:            |  done
*/
int meta_close(){

    /*global*/ meta_end = clock();
    double differential_time = (/*global*/ meta_end - /*global*/ meta_start) / CLOCKS_PER_SEC;
    char endstring[102];

    if(meta_write_int_0D("FINISHED_WALLTIME", NULL, differential_time))
        return 1;
    
    for(size_t i = 0; i < 101; i++)
        endstring[i] = '-';
    endstring[101] = '\0';

    fprintf(/*global*/ fh_meta_out, "%s\n", endstring);

    if(/*global*/ fh_meta_in != NULL)
        fclose(/*global*/ fh_meta_in);
    if(/*global*/ fh_meta_out != NULL)
        fclose(/*global*/ fh_meta_out);

    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method to get the number of lines in a file.
*
* \param[in] fh File handle to operate on.
* \return -1 if the file contains no newlines, else the number of lines.
*/
size_t meta_count_lines(FILE *fh){
    if(fh == NULL) 
        return -1;

    char ch;
    bool touched = false;
    size_t number_of_lines = 0, i;
    
    for(i = 0; (ch = fgetc(fh)) != EOF; i++)
        if(ch == '\n'){
            number_of_lines++;
            touched = true;
        }
    if(fseek(fh, 1 - i, SEEK_CUR)) 
        return -1;
    return touched ? number_of_lines : -1;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method to parse the basename.
*
* \param[in] filename Full name of the file.
* \param[in] suf Expected suffix.
*/
int meta_parse_basename(char *filename, char *suf){
    if(filename == NULL || suf == NULL)
        return 1;

    size_t number_of_lines, counter = 0, errortext_length = strlen(/*global*/ in.full_name) + 72;
    char *filename_part, *basename_end, buffer[META_MCL], *ptr;
    char errortext[errortext_length], full_name_buffer[META_MCL];
    bool error = false;
    int err_number;

    if(strnlen(filename, META_MCL) == META_MCL) 
        return 1;
    strcpy(full_name_buffer, filename);
    __meta_zero_basename_struct((void *) /*global*/ &in);
    strcpy(/*global*/ in.full_name, full_name_buffer);
    
    filename_part = strtok(full_name_buffer, "/");
    if(filename_part == NULL) 
        return 1;
    
    ptr = /*global*/ in.path;
    while(filename_part != NULL){
        if(strstr(filename_part, suf) != NULL)
            break;
        else{
            ptr = __meta_fast_strcat(ptr, "/");
            ptr = __meta_fast_strcat(ptr, filename_part);
            if(ptr == NULL)
                return 1;
            filename_part = strtok(NULL, "/");
        }
    }
    if(filename_part == NULL) 
        return 1;

    strcpy(/*global*/ in.basename, filename_part);
    basename_end = strstr(/*global*/ in.basename, suf);
    *basename_end = '\0';
    strcat(/*global*/ in.path, "/");
    strcpy(/*global*/ in.path_and_basename, /*global*/ in.path);
    strcat(/*global*/ in.path_and_basename, /*global*/ in.basename);
    
    strcpy(buffer, /*global*/ in.basename);
    filename_part = strtok(buffer, META_BASENAME_SEPARATOR);
    if(filename_part == NULL){
        err_number = 0;
        error = true;
    }
    if(!error)
        strcpy(/*global*/ in.dataset, filename_part);
    filename_part = strtok(NULL, META_BASENAME_SEPARATOR);
    if(filename_part == NULL){
        if(!error)
            err_number = 1;
        error = true;
    }
    if(!error)
        strcpy(/*global*/ in.type, filename_part);
    filename_part = strtok(NULL, META_BASENAME_SEPARATOR);
    if(filename_part == NULL){
        if(!error)
            err_number = 2;
        error = true;
    }
    if(!error)
        strcpy(/*global*/ in.purpose, filename_part);
    filename_part = strtok(NULL, META_BASENAME_SEPARATOR);
    if(filename_part == NULL){
        if(!error)
            err_number = 3;
        error = true;
    }
    if(!error)
        strcpy(/*global*/ in.app, filename_part);
    filename_part = strtok(NULL, META_BASENAME_SEPARATOR);
    if(filename_part == NULL){
        if(!error)
            err_number = 4;
        error = true;
    }
    if(!error)
        strcpy(/*global*/ in.features, filename_part);
    
    filename_part = strtok(NULL, META_BASENAME_SEPARATOR);
    
    if(filename_part != NULL || error){
        if(error)
            counter = err_number;
        else if(filename_part != NULL){
            counter = 5;
            while(filename_part != NULL){
                counter++;
                filename_part = strtok(NULL, META_BASENAME_SEPARATOR);
            }
        }
        snprintf(
            errortext, errortext_length, 
            "Basename with %zu segments given, which is invalid: %s", 
            counter, /*global*/ in.path_and_basename
        );
        return __meta_print_error(stdout, errortext, 1);
    }
        
    
    strcpy(/*global*/ out.full_name, /*global*/ in.full_name);
    strcpy(/*global*/ out.path, /*global*/ in.path);
    strcpy(/*global*/ out.path_and_basename, /*global*/ in.path_and_basename);
    strcpy(/*global*/ out.basename, /*global*/ in.basename);
    strcpy(/*global*/ out.dataset, /*global*/ in.dataset);
    strcpy(/*global*/ out.type, /*global*/ in.type);
    strcpy(/*global*/ out.purpose, /*global*/ in.purpose);
    strcpy(/*global*/ out.app, /*global*/ in.app);
    strcpy(/*global*/ out.features, /*global*/ in.features);
    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method to check and in case trunctate length of given unit.
*
* \param[in] unit The physical unit to be checked.
* \return 0 if ok, 1 if trunctated, on error -1.
*/
/*
* UNIFY HEAD:             |  fortran error (name_fortran=check_unit)
* UNIFY DOC:              |  fortran error (unprecise,wrong)
* FUNCTIONAL (THEORY):    |  done, fortran error(implementation_error=true["does not trunctate the keyword"])
* FUNCTIONAL (PRACTICE):  |  
* UNIT TESTING:           | 
* ADD  ERROR HANDLER:     |  done
* ADD  INPUT SANITIZER:   |  done
* MAKE PRETTY:            |  done
*/
int meta_check_unit(char *unit){
    if(unit == NULL) 
        return -1;
    if(strnlen(unit, META_UCL) < META_UCL) 
        return 0;
    
    char errortext[strlen(unit) + 74], *ptr;
    ptr = __meta_fast_strcat(errortext, "The unit ");
    ptr = __meta_fast_strcat(ptr, unit);
    ptr = __meta_fast_strcat(ptr, "is longer than\nthe convention allows and therefore truncated!\n");

    fprintf(stdout, errortext);

    unit[META_KCL - 1] = '\0';
    return 1;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method to check and in case trunctate length of given keyword.
*
* \param[in] keyword The keyword to be checked.
* \return 0 if ok, 1 if trunctated, on error -1.
*/
/*
* UNIFY HEAD:             |  fortran error (name_fortran=check_keyword)
* UNIFY DOC:              |  fortran error (unprecise)
* FUNCTIONAL (THEORY):    |  done
* FUNCTIONAL (PRACTICE):  | 
* UNIT TESTING:           | 
* ADD  ERROR HANDLER:     |  done
* ADD  INPUT SANITIZER:   |  done
* MAKE PRETTY:            |  done
*/
int meta_check_keyword(char *keyword){
    if(keyword == NULL) 
        return -1;
    if(strnlen(keyword, META_KCL) < META_KCL) 
        return 0;
    
    char errortext[strlen(keyword) + 77], *ptr;
    ptr = __meta_fast_strcat(errortext, "The keyword ");
    ptr = __meta_fast_strcat(ptr, keyword);
    ptr = __meta_fast_strcat(ptr, "is longer than\nthe convention allows and therefore truncated!\n");

    fprintf(stdout, errortext);

    keyword[META_KCL - 1] = '\0';
    return 1;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method to take the sha256 hashing function of a certain binary and write it to metafile.
*
* \param[in] binary_name The file to be hashed.
* \return 0 on success, 1 otherwise.
*/
/*
* UNIFY HEAD:             |  done
* UNIFY DOC:              |  fortran error
* FUNCTIONAL (THEORY):    |  done
* FUNCTIONAL (PRACTICE):  | 
* UNIT TESTING:           | 
* ADD  ERROR HANDLER:     | 
* ADD  INPUT SANITIZER:   |  done
* MAKE PRETTY:            |  done
*/
int meta_write_sha256sum(char *binary_name){
    if(binary_name == NULL) 
        return 1;

    char keyword[] = "SHA256SUM_OF_BINARY";
    size_t hash_size = 258 * sizeof(char);
    size_t command_size = strlen(binary_name + 46) * sizeof(char);
    char hash[hash_size], command[command_size], *ptr;
    int error = 0;
    FILE *fh;

    if(__meta_get_filesize("temp_buffer") != -1)
        if(system("rm -r temp_buffer")) 
            return 1;

    if(
        system("which cut > /dev/null 2> /dev/null") ||
        system("which sha256sum > /dev/null 2> /dev/null")
    ) 
        return 1;
    
    fh = fopen("temp_buffer", "w");
    if(fh == NULL) 
        return 1;
    fclose(fh);
    
    ptr = __meta_fast_strcat(command, "sha256sum ");
    ptr = __meta_fast_strcat(ptr, binary_name);
    ptr = __meta_fast_strcat(ptr, " | cut -d ' ' -f 1 >> 'temp_buffer'");
    if(ptr == NULL)
        return 1;

    if(!system(command)){
        fh = fopen("temp_buffer", "r");
        if(fh == NULL) 
            return 1;
        if(fgets(hash, hash_size, fh) == NULL) 
            error = 1;
        fclose(fh);
        
        if(!error)
            if(meta_write_string(keyword, hash)) 
                error = 1;
    }
    
    if(__meta_get_filesize("temp_buffer") != -1)
        if(system("rm -r temp_buffer")) 
            return 1;

    return error;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief Method to check filesize and delete file if 0.
*
* \param[in] filename Path and name of the file to be checked and deleted.
* \return 0 if file was deleted, 1 if not. On error -1.
*/
/*
* UNIFY HEAD:             |  done
* UNIFY DOC:              |  done
* FUNCTIONAL (THEORY):    |  done
* FUNCTIONAL (PRACTICE):  | 
* UNIT TESTING:           | 
* ADD  ERROR HANDLER:     | 
* ADD  INPUT SANITIZER:   |  done
* MAKE PRETTY:            |  done
*/
int meta_delete_empty_file(char *filename){
    if(filename == NULL) 
        return -1;

    ssize_t filesize = __meta_get_filesize(filename);
    char command[strlen(filename) + 7], *ptr;

    if(filesize == -1) 
        return -1;
    if(filesize > 0) 
        return 1;

    ptr = __meta_fast_strcat(command, "rm -r ");
    ptr = __meta_fast_strcat(ptr, filename);
    if(ptr == NULL)
        return -1;

    if(system(command)) 
        return -1;
    return 0;
}

// ---- private function implementations ----

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief library setter to update a line in a metafile structure.
*
* \private
*
* \param[in] metafile A pointer to a metafile structure.
* \param[in] line_number The line one wants to alter.
* \param[in] string The replacement string for the requested line.
* \return on success 0, 1 otherwise.
*/
static int __meta_set_metafile_string(metafile *metafile, unsigned int line_number, char *string){
    if(metafile == NULL || string == NULL) 
        return 1;
    if(line_number >= metafile -> number_of_lines) 
        return 1;
    if(strnlen(string, META_MCL) == META_MCL) 
        return 1;

    strcpy(metafile -> content + line_number * META_MCL, string);
    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief library getter to obtain a reference (raw pointer) to a metafile line (= a string).
*
* \private
*
* \param[in] metafile A pointer to a metafile structure.
* \param[in] line_number The line one wants to get a reference on.
* \return Pointer to the requested line if possible, else or on error NULL.
*/
static char *__meta_get_metafile_string_reference(metafile *metafile, unsigned int line_number){
    if(metafile == NULL) 
        return NULL;
    if(line_number >= metafile -> number_of_lines) 
        return NULL;
    return (metafile -> content) + line_number * META_MCL * sizeof(char);
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief library getter to obtain a safe copy of a metafile line (= a string).
*
* \private
*
* \param[in] metafile A pointer to a metafile structure.
* \param[in] line_number The line one wants to get a copy of.
* \param[out] copy The buffer the metafile line should be copied to.
* \return on success 0, else 1.
*/
static int __meta_get_metafile_string_copy(metafile *metafile, unsigned int line_number, char *copy){
    if(metafile == NULL || copy == NULL) 
        return 1;
    if(line_number >= metafile -> number_of_lines) 
        return 1;
    if(strnlen(metafile -> content + line_number * META_MCL, META_MCL) == META_MCL) 
        return 1;

    strcpy(copy, metafile -> content + line_number * META_MCL);
    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief small library iternal to set all elements in a basename structure to zero
*
* \private
*
* \param[out] basename A pointer to the basename struct.
*/
static void __meta_zero_basename_struct(basename *basename){
        basename -> full_name[0] = '\0';
        basename -> path_and_basename[0] = '\0';
        basename -> basename[0] = '\0';
        basename -> dataset[0] = '\0';
        basename -> type[0] = '\0';
        basename -> purpose[0] = '\0';
        basename -> app[0] = '\0';
        basename -> features[0] = '\0';
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief small library internal to zero all elements in an array.
*
* \private
*
* \param[in] array_size The size of the array in bytes.
* \param[inout] array The array to zero.
*/
static void __meta_zero_array(size_t array_size, void *array){
    if(array == NULL) 
        return;
    char *internal_array = (char *) array; 
    for(int i = 0; i < array_size; i++) 
        internal_array[i] = '\0';
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief small library internal to get the filesize and existence of a file.
* 
* \private
*
* \param[in] filename The filename to be checked.
* \return -1 if file does not exist, filesize in bytes otherwise.
*/
static ssize_t __meta_get_filesize(char *filename){
    if(filename == NULL) 
        return -1;
    struct stat st;
    if(stat(filename, &st) != 0) 
        return -1; 
    return (ssize_t) st.st_size;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief library internal to get the file extension (suffix).
*
* \private
*
* \param[in] filename The filename to be parsed.
* \param[out] buffer The buffer to store the filename extension without META_SUFFIX_DECLARATOR.
* \return 1 on error, on success 0.
*/
static int __meta_get_file_suffix(char *filename_with_suffix, char *buffer){
    if(filename_with_suffix == NULL) 
        return 1;
    if(buffer == NULL)
        return 1;

    size_t filename_length = strlen(filename_with_suffix);
    char *old_token, *token, internal_buffer[filename_length + 1];

    strcpy(internal_buffer, filename_with_suffix);
    
    token = strtok(internal_buffer, META_SUFFIX_DECLARATOR);
    if(token == NULL) 
        return 1;
    while(token != NULL){
        old_token = token;
        token = strtok(NULL, META_SUFFIX_DECLARATOR);
    }
    strcpy(buffer, old_token);
    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief library internal to concatenate multiple strings in an efficient way.
*
* Description:
*   Always use large enough destination buffer!
*   Usage is: 
*   __meta_fast_strcat(... (__meta_fast_strcat(destination, source1), ...), sourceN)
*   Trick is: this implementation does not search over the string for \0, but is
*   more a mixture of strcat and strcpy, for it is used as strcat, but implemented
*   similar to strcpy.
*
* \private
*
* \param[out] destination The destination buffer to concatenate to. Has to be large enough.
* \param[in] source. The source buffer from where to concatenate.
* \return NULL on error, on success pointer to the current end ('\0') of destination.
*/
static char* __meta_fast_strcat(char *destination, char *source){
    if(destination == NULL || source == NULL) 
        return NULL;

    while(*source != '\0')
        *destination++ = *source++;
    *destination = '\0';
    return destination;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief library internal to replace each occurrence of a char in a string with another.
*
* Description:
*   If replacement == '\0', the method will terminate after the first occurrence and
*   replacement of to_be_replaced, as it is unnessesary to execute replacement further
*   (The string is already terminated).
*
* \private
*
* \param[inout] string The string to be altered.
* \param[in] to_be_replaced Char to be replaced in string.
* \param[in] replacement. Replacement char.
*/
static int __meta_replace_char(char *string, char to_be_replaced, char replacement){
    if(string == NULL)
        return 1;

    for(size_t i = 0; string[i] != '\0'; i++)
        if(string[i] == to_be_replaced){
            string[i] = replacement;
            if(replacement == '\0')
                break;
        }
    return 0;
}

/**
* \author Jonathan Schäfer, hpcjscha@hlrs.de, HLRS/NUM
*
* \brief library internal to print an error to a certain file.
*
* Description:
*   Returns the given error value. Can be used in this way:
*   return __meta_print_error(fh, "Error occurred", 1); .
*
* \private
*
* \param[in] fh The file handle to write the error to.
* \param[in] message The error message. Nullable for generic message.
* \param[in] error Integer error value.
* \return The same integer given with "error".
*/
static int __meta_print_error(FILE *fh, char *message, int error){
    if(fh == NULL)
        return 1;

    char errormessage[message != NULL ? strlen(message) + 1 : 61];
    if(message != NULL)
        strcpy(errormessage, message);
    else
        strcpy(errormessage, "Generic error occurred during metafile execution!");
    
    fprintf(fh, "%s WITH ERRORCODE %d\n", errormessage, error);
    return error;
}

