/*
 * Copyright 1993-2002 Christopher Seiwald and Perforce Software, Inc.
 *
 * This file is part of Jam - see jam.c for Copyright information.
 */

#include "jam.h"

#include "lists.h"
#include "parse.h"
#include "builtins.h"
#include "rules.h"
#include "filesys.h"
#include "object.h"
#include "regexp.h"
#include "frames.h"
#include "hash.h"
#include "strings.h"
#include "pwd.h"
#include "pathsys.h"
#include "make.h"
#include "hdrmacro.h"
#include "compile.h"
#include "native.h"
#include "variable.h"
#include "timestamp.h"
#include "md5.h"
#include "constants.h"
#include <ctype.h>

#if defined(USE_EXECUNIX)
# include <sys/types.h>
# include <sys/wait.h>
#else
/*
  NT does not have wait() and associated macros, it uses the return value
  of system() instead. Status code group are documented at
  http://msdn.microsoft.com/en-gb/library/ff565436.aspx
*/
# define WIFEXITED(w)  (((w) & 0XFFFFFF00) == 0)
# define WEXITSTATUS(w)(w)
#endif

/*
 * builtins.c - builtin jam rules
 *
 * External routines:
 *
 *  load_builtin() - define builtin rules
 *
 * Internal routines:
 *
 *  builtin_depends() - DEPENDS/INCLUDES rule.
 *  builtin_echo() - ECHO rule.
 *  builtin_exit() - EXIT rule.
 *  builtin_flags() - NOCARE, NOTFILE, TEMPORARY rule.
 *  builtin_glob() - GLOB rule.
 *  builtin_match() - MATCH rule.
 *
 * 01/10/01 (seiwald) - split from compile.c
 */


/*
 * compile_builtin() - define builtin rules
 */

#define P0 (PARSE *)0
#define C0 (OBJECT *)0

#if defined( OS_NT ) || defined( OS_CYGWIN )
    LIST * builtin_system_registry      ( FRAME *, int );
    LIST * builtin_system_registry_names( FRAME *, int );
#endif

int glob( const char * s, const char * c );

void backtrace        ( FRAME * );
void backtrace_line   ( FRAME * );
void print_source_line( FRAME * );


RULE * bind_builtin( const char * name_, LIST * (* f)( FRAME *, int flags ), int flags, const char * * args )
{
    FUNCTION * func;
    RULE * result;
    argument_list* arg_list = 0;
    OBJECT * name = object_new( name_ );

    if ( args )
    {
        arg_list = args_new();
        lol_build( arg_list->data, args );
    }

    func = function_builtin( f, flags );

    result = new_rule_body( root_module(), name, arg_list, func, 1 );

    function_free( func );

    object_free( name );

    return result;
}


RULE * duplicate_rule( const char * name_, RULE * other )
{
    OBJECT * name = object_new( name_ );
    RULE * result = import_rule( other, root_module(), name );
    object_free( name );
    return result;
}


void load_builtins()
{
    duplicate_rule( "Always",
      bind_builtin( "ALWAYS",
                    builtin_flags, T_FLAG_TOUCHED, 0 ) );

    duplicate_rule( "Depends",
      bind_builtin( "DEPENDS",
                    builtin_depends, 0, 0 ) );

    duplicate_rule( "echo",
    duplicate_rule( "Echo",
      bind_builtin( "ECHO",
                    builtin_echo, 0, 0 ) ) );

    {
        const char * args[] = { "message", "*", ":", "result-value", "?", 0 };
        duplicate_rule( "exit",
        duplicate_rule( "Exit",
          bind_builtin( "EXIT",
                        builtin_exit, 0, args ) ) );
    }

    {
        const char * args[] = { "directories", "*", ":", "patterns", "*", ":", "case-insensitive", "?", 0 };
        duplicate_rule( "Glob",
                        bind_builtin( "GLOB", builtin_glob, 0, args ) );
    }

    {
        const char * args[] = { "patterns", "*", 0 };
        bind_builtin( "GLOB-RECURSIVELY",
                      builtin_glob_recursive, 0, args );
    }

    duplicate_rule( "Includes",
      bind_builtin( "INCLUDES",
                    builtin_depends, 1, 0 ) );

    {
        const char * args[] = { "targets", "*", ":", "targets-to-rebuild", "*", 0 };
        bind_builtin( "REBUILDS",
                      builtin_rebuilds, 0, args );
    }

    duplicate_rule( "Leaves",
      bind_builtin( "LEAVES",
                    builtin_flags, T_FLAG_LEAVES, 0 ) );

    duplicate_rule( "Match",
      bind_builtin( "MATCH",
                    builtin_match, 0, 0 ) );

    {
        const char * args[] = { "string", ":", "delimiters" };
        bind_builtin( "SPLIT_BY_CHARACTERS", 
                      builtin_split_by_characters, 0, 0 );
    }

    duplicate_rule( "NoCare",
      bind_builtin( "NOCARE",
                    builtin_flags, T_FLAG_NOCARE, 0 ) );

    duplicate_rule( "NOTIME",
    duplicate_rule( "NotFile",
      bind_builtin( "NOTFILE",
                    builtin_flags, T_FLAG_NOTFILE, 0 ) ) );

    duplicate_rule( "NoUpdate",
      bind_builtin( "NOUPDATE",
                    builtin_flags, T_FLAG_NOUPDATE, 0 ) );

    duplicate_rule( "Temporary",
      bind_builtin( "TEMPORARY",
                    builtin_flags, T_FLAG_TEMP, 0 ) );

      bind_builtin( "ISFILE",
                    builtin_flags, T_FLAG_ISFILE, 0 );

    duplicate_rule( "HdrMacro",
      bind_builtin( "HDRMACRO",
                    builtin_hdrmacro, 0, 0 ) );

    /* FAIL_EXPECTED is used to indicate that the result of a target build
     * action should be inverted (ok <=> fail) this can be useful when
     * performing test runs from Jamfiles.
     */
      bind_builtin( "FAIL_EXPECTED",
                    builtin_flags, T_FLAG_FAIL_EXPECTED, 0 );

      bind_builtin( "RMOLD",
                    builtin_flags, T_FLAG_RMOLD, 0 );

      {
          const char * args[] = { "targets", "*", 0 };
          bind_builtin( "UPDATE",
                        builtin_update, 0, args );
      }

      {
          const char * args[] = { "targets", "*", 
                            ":", "log", "?",
                            ":", "ignore-minus-n", "?", 
                            ":", "ignore-minus-q", "?", 0 };
          bind_builtin( "UPDATE_NOW",
                        builtin_update_now, 0, args );
      }

      {
          const char * args[] = { "string", "pattern", "replacements", "+", 0 };
          duplicate_rule( "subst",
            bind_builtin( "SUBST",
                          builtin_subst, 0, args ) );
      }

      {
          const char * args[] = { "module", "?", 0 };
          bind_builtin( "RULENAMES",
                         builtin_rulenames, 0, args );
      }


      {
          const char * args[] = { "module", "?", 0 };
          bind_builtin( "VARNAMES",
                         builtin_varnames, 0, args );
      }

      {
          const char * args[] = { "module", "?", 0 };
          bind_builtin( "DELETE_MODULE",
                         builtin_delete_module, 0, args );
      }

      {
          const char * args[] = { "source_module", "?",
                            ":", "source_rules", "*",
                            ":", "target_module", "?",
                            ":", "target_rules", "*",
                            ":", "localize", "?", 0 };
          bind_builtin( "IMPORT",
                        builtin_import, 0, args );
      }

      {
          const char * args[] = { "module", "?", ":", "rules", "*", 0 };
          bind_builtin( "EXPORT",
                        builtin_export, 0, args );
      }

      {
          const char * args[] = { "levels", "?", 0 };
          bind_builtin( "CALLER_MODULE",
                         builtin_caller_module, 0, args );
      }

      {
          const char * args[] = { "levels", "?", 0 };
          bind_builtin( "BACKTRACE",
                        builtin_backtrace, 0, args );
      }

      {
          const char * args[] = { 0 };
          bind_builtin( "PWD",
                        builtin_pwd, 0, args );
      }

      {
          const char * args[] = { "target", "*", ":", "path", "*", 0 };
          bind_builtin( "SEARCH_FOR_TARGET",
                        builtin_search_for_target, 0, args );
      }

      {
          const char * args[] = { "modules_to_import", "+", ":", "target_module", "?", 0 };
          bind_builtin( "IMPORT_MODULE",
                        builtin_import_module, 0, args );
      }

      {
          const char * args[] = { "module", "?", 0 };
          bind_builtin( "IMPORTED_MODULES",
                        builtin_imported_modules, 0, args );
      }

      {
          const char * args[] = { "instance_module", ":", "class_module", 0 };
          bind_builtin( "INSTANCE",
                        builtin_instance, 0, args );
      }

      {
          const char * args[] = { "sequence", "*", 0 };
          bind_builtin( "SORT",
                        builtin_sort, 0, args );
      }

      {
          const char * args[] = { "path_parts", "*", 0 };
          bind_builtin( "NORMALIZE_PATH",
                        builtin_normalize_path, 0, args );
      }

      {
          const char * args[] = { "args", "*", 0 };
          bind_builtin( "CALC",
                        builtin_calc, 0, args );
      }

      {
          const char * args[] = { "module", ":", "rule", 0 };
          bind_builtin( "NATIVE_RULE",
                        builtin_native_rule, 0, args );
      }

      {
          const char * args[] = { "module", ":", "rule", ":", "version", 0 };
          bind_builtin( "HAS_NATIVE_RULE",
                        builtin_has_native_rule, 0, args );
      }

      {
          const char * args[] = { "module", "*", 0 };
          bind_builtin( "USER_MODULE",
                        builtin_user_module, 0, args );
      }

      {
          const char * args[] = { 0 };
          bind_builtin( "NEAREST_USER_LOCATION",
                        builtin_nearest_user_location, 0, args );
      }

      {
          const char * args[] = { "file", 0 };
          bind_builtin( "CHECK_IF_FILE",
                        builtin_check_if_file, 0, args );
      }

#ifdef HAVE_PYTHON
      {
          const char * args[] = { "python-module", ":", "function", ":",
                            "jam-module", ":", "rule-name", 0 };
          bind_builtin( "PYTHON_IMPORT_RULE",
                        builtin_python_import_rule, 0, args );
      }
#endif

# if defined( OS_NT ) || defined( OS_CYGWIN )
      {
          const char * args[] = { "key_path", ":", "data", "?", 0 };
          bind_builtin( "W32_GETREG",
                        builtin_system_registry, 0, args );
      }

      {
          const char * args[] = { "key_path", ":", "result-type", 0 };
          bind_builtin( "W32_GETREGNAMES",
                        builtin_system_registry_names, 0, args );
      }
# endif

      {
          const char * args[] = { "command", ":", "*", 0 };
          duplicate_rule( "SHELL",
            bind_builtin( "COMMAND",
                          builtin_shell, 0, args ) );
      }

      {
          const char * args[] = { "string", 0 };
          bind_builtin( "MD5",
                        builtin_md5, 0, args ) ;
      }

      {
          const char * args[] = { "name", ":", "mode", 0 };
          bind_builtin( "FILE_OPEN",
                        builtin_file_open, 0, args );
      }

      {
          const char * args[] = { "string", ":", "width", 0 };
          bind_builtin( "PAD",
                        builtin_pad, 0, args );
      }

      {
          const char * args[] = { "targets", "*", 0 };
          bind_builtin( "PRECIOUS",
                        builtin_precious, 0, args );
      }

      {
          const char * args [] = { 0 };
          bind_builtin( "SELF_PATH", builtin_self_path, 0, args );
      }

      {
          const char * args [] = { "path", 0 };
          bind_builtin( "MAKEDIR", builtin_makedir, 0, args );
      }

      {
          const char * args[] = { "directories", "*", 0 };
          bind_builtin( "RESCAN",
                        builtin_rescan, 0, args );
      }

      /* Initialize builtin modules. */
      init_set();
      init_path();
      init_regex();
      init_property_set();
      init_sequence();
      init_order();
      init_svnrev();
}


/*
 * builtin_calc() - CALC rule.
 *
 * The CALC rule performs simple mathematical operations on two arguments.
 */

LIST * builtin_calc( FRAME * frame, int flags )
{
    LIST * arg = lol_get( frame->args, 0 );

    LIST * result = 0;
    long lhs_value;
    long rhs_value;
    long result_value;
    char buffer [ 16 ];
    char const * lhs;
    char const * op;
    char const * rhs;

    if ( arg == 0 ) return L0;
    lhs = object_str( arg->value );

    arg = list_next( arg );
    if ( arg == 0 ) return L0;
    op = object_str( arg->value );

    arg = list_next( arg );
    if ( arg == 0 ) return L0;
    rhs = object_str( arg->value );

    lhs_value = atoi( lhs );
    rhs_value = atoi( rhs );

    if ( strcmp( "+", op ) == 0 )
    {
        result_value = lhs_value + rhs_value;
    }
    else if ( strcmp( "-", op ) == 0 )
    {
        result_value = lhs_value - rhs_value;
    }
    else
    {
        return L0;
    }

    sprintf( buffer, "%ld", result_value );
    result = list_new( result, object_new( buffer ) );
    return result;
}


/*
 * builtin_depends() - DEPENDS/INCLUDES rule.
 *
 * The DEPENDS/INCLUDES builtin rule appends each of the listed sources on the
 * dependency/includes list of each of the listed targets. It binds both the
 * targets and sources as TARGETs.
 */

LIST * builtin_depends( FRAME * frame, int flags )
{
    LIST * targets = lol_get( frame->args, 0 );
    LIST * sources = lol_get( frame->args, 1 );
    LIST * l;

    for ( l = targets; l; l = list_next( l ) )
    {
        TARGET * t = bindtarget( l->value );

        /* If doing INCLUDES, switch to the TARGET's include */
        /* TARGET, creating it if needed.  The internal include */
        /* TARGET shares the name of its parent. */

        if ( flags )
        {
            if ( !t->includes )
            {
                t->includes = copytarget( t );
                t->includes->original_target = t;
            }
            t = t->includes;
        }

        t->depends = targetlist( t->depends, sources );
    }

    /* Enter reverse links */
    for ( l = sources; l; l = list_next( l ) )
    {
        TARGET * s = bindtarget( l->value );
        s->dependants = targetlist( s->dependants, targets );
    }

    return L0;
}


/*
 * builtin_rebuilds() - REBUILDS rule.
 *
 * The REBUILDS builtin rule appends each of the listed rebuild-targets in its
 * 2nd argument on the rebuilds list of each of the listed targets in its first
 * argument.
 */

LIST * builtin_rebuilds( FRAME * frame, int flags )
{
    LIST * targets = lol_get( frame->args, 0 );
    LIST * rebuilds = lol_get( frame->args, 1 );
    LIST * l;

    for ( l = targets; l; l = list_next( l ) )
    {
        TARGET * t = bindtarget( l->value );
        t->rebuilds = targetlist( t->rebuilds, rebuilds );
    }

    return L0;
}


/*
 * builtin_echo() - ECHO rule.
 *
 * The ECHO builtin rule echoes the targets to the user. No other actions are
 * taken.
 */

LIST * builtin_echo( FRAME * frame, int flags )
{
    list_print( lol_get( frame->args, 0 ) );
    printf( "\n" );
    fflush( stdout );
    return L0;
}


/*
 * builtin_exit() - EXIT rule.
 *
 * The EXIT builtin rule echoes the targets to the user and exits the program
 * with a failure status.
 */

LIST * builtin_exit( FRAME * frame, int flags )
{
    list_print( lol_get( frame->args, 0 ) );
    printf( "\n" );
    if ( lol_get( frame->args, 1 ) )
    {
        exit( atoi( object_str( lol_get( frame->args, 1 )->value ) ) );
    }
    else
    {
        exit( EXITBAD );  /* yeech */
    }
    return L0;
}


/*
 * builtin_flags() - NOCARE, NOTFILE, TEMPORARY rule.
 *
 * Builtin_flags() marks the target with the appropriate flag, for use by make0().
 * It binds each target as a TARGET.
 */

LIST * builtin_flags( FRAME * frame, int flags )
{
    LIST * l = lol_get( frame->args, 0 );
    for ( ; l; l = list_next( l ) )
        bindtarget( l->value )->flags |= flags;
    return L0;
}


/*
 * builtin_globbing() - GLOB rule.
 */

struct globbing
{
    LIST * patterns;
    LIST * results;
    LIST * case_insensitive;
};


static void downcase_inplace( char * p )
{
    for ( ; *p; ++p )
        *p = tolower( *p );
}


static void builtin_glob_back
(
    void   * closure,
    OBJECT * file,
    int      status,
    time_t   time
)
{
    PROFILE_ENTER( BUILTIN_GLOB_BACK );

    struct globbing * globbing = (struct globbing *)closure;
    LIST            * l;
    PATHNAME          f;
    string            buf[ 1 ];

    /* Null out directory for matching. We wish we had file_dirscan() pass up a
     * PATHNAME.
     */
    path_parse( object_str( file ), &f );
    f.f_dir.len = 0;

    /* For globbing, we unconditionally ignore current and parent directory
     * items. Since these items always exist, there is no reason why caller of
     * GLOB would want to see them. We could also change file_dirscan(), but
     * then paths with embedded "." and ".." would not work anywhere.
    */
    if ( !strcmp( f.f_base.ptr, "." ) || !strcmp( f.f_base.ptr, ".." ) )
    {
        PROFILE_EXIT( BUILTIN_GLOB_BACK );
        return;
    }

    string_new( buf );
    path_build( &f, buf, 0 );

    if ( globbing->case_insensitive )
        downcase_inplace( buf->value );

    for ( l = globbing->patterns; l; l = l->next )
    {
        if ( !glob( object_str( l->value ), buf->value ) )
        {
            globbing->results = list_new( globbing->results, object_copy( file ) );
            break;
        }
    }

    string_free( buf );

    PROFILE_EXIT( BUILTIN_GLOB_BACK );
}


static LIST * downcase_list( LIST * in )
{
    LIST * result = 0;

    string s[ 1 ];
    string_new( s );

    while ( in )
    {
        string_copy( s, object_str( in->value ) );
        downcase_inplace( s->value );
        result = list_append( result, list_new( 0, object_new( s->value ) ) );
        in = in->next;
    }

    string_free( s );
    return result;
}


LIST * builtin_glob( FRAME * frame, int flags )
{
    LIST * l = lol_get( frame->args, 0 );
    LIST * r = lol_get( frame->args, 1 );

    struct globbing globbing;

    globbing.results = L0;
    globbing.patterns = r;

    globbing.case_insensitive
# if defined( OS_NT ) || defined( OS_CYGWIN )
       = l;  /* Always case-insensitive if any files can be found. */
# else
       = lol_get( frame->args, 2 );
# endif

    if ( globbing.case_insensitive )
        globbing.patterns = downcase_list( r );

    for ( ; l; l = list_next( l ) )
        file_dirscan( l->value, builtin_glob_back, &globbing );

    if ( globbing.case_insensitive )
        list_free( globbing.patterns );

    return globbing.results;
}


static int has_wildcards( char const * str )
{
    size_t const index = strcspn( str, "[]*?" );
    return str[ index ] == '\0' ? 0 : 1;
}


/*
 * If 'file' exists, append 'file' to 'list'. Returns 'list'.
 */

static LIST * append_if_exists( LIST * list, OBJECT * file )
{
    time_t time;
    timestamp( file, &time );
    return time > 0
        ? list_new( list, object_copy( file ) )
        : list;
}


LIST * glob1( OBJECT * dirname, OBJECT * pattern )
{
    LIST * plist = list_new( L0, object_copy(pattern) );
    struct globbing globbing;

    globbing.results = L0;
    globbing.patterns = plist;

    globbing.case_insensitive
# if defined( OS_NT ) || defined( OS_CYGWIN )
       = plist;  /* always case-insensitive if any files can be found */
# else
       = L0;
# endif

    if ( globbing.case_insensitive )
        globbing.patterns = downcase_list( plist );

    file_dirscan( dirname, builtin_glob_back, &globbing );

    if ( globbing.case_insensitive )
        list_free( globbing.patterns );

    list_free( plist );

    return globbing.results;
}


LIST * glob_recursive( const char * pattern )
{
    LIST * result = L0;

    /* Check if there's metacharacters in pattern */
    if ( !has_wildcards( pattern ) )
    {
        /* No metacharacters. Check if the path exists. */
        OBJECT * p = object_new( pattern );
        result = append_if_exists( result, p );
        object_free( p );
    }
    else
    {
        /* Have metacharacters in the pattern. Split into dir/name. */
        PATHNAME path[ 1 ];
        path_parse( pattern, path );

        if ( path->f_dir.ptr )
        {
            LIST * dirs = L0;
            string dirname[ 1 ];
            string basename[ 1 ];
            string_new( dirname );
            string_new( basename );

            string_append_range( dirname, path->f_dir.ptr,
                                path->f_dir.ptr + path->f_dir.len );

            path->f_grist.ptr = 0;
            path->f_grist.len = 0;
            path->f_dir.ptr = 0;
            path->f_dir.len = 0;
            path_build( path, basename, 0 );

            dirs =  has_wildcards( dirname->value )
                ? glob_recursive( dirname->value )
                : list_new( dirs, object_new( dirname->value ) );

            if ( has_wildcards( basename->value ) )
            {
                LIST * d;
                OBJECT * b = object_new( basename->value );
                for ( d = dirs ; d; d = d->next )
                    result = list_append( result, glob1( d->value, b ) );
                object_free( b );
            }
            else
            {
                LIST * d;
                string file_string[ 1 ];
                string_new( file_string );

                /* No wildcard in basename. */
                for ( d = dirs ; d; d = d->next )
                {
                    OBJECT * p;
                    path->f_dir.ptr = object_str( d->value );
                    path->f_dir.len = strlen( object_str( d->value ) );
                    path_build( path, file_string, 0 );

                    p = object_new( file_string->value );

                    result = append_if_exists( result, p );

                    object_free( p );

                    string_truncate( file_string, 0 );
                }

                string_free( file_string );
            }

            string_free( dirname );
            string_free( basename );

            list_free( dirs );
        }
        else
        {
            /** No directory, just a pattern. */
            OBJECT * d = object_new( "." );
            OBJECT * p = object_new( pattern );
            result = list_append( result, glob1( d, p ) );
            object_free( p );
            object_free( d );
        }
    }

    return result;
}


LIST * builtin_glob_recursive( FRAME * frame, int flags )
{
    LIST * result = L0;
    LIST * l = lol_get( frame->args, 0 );
    for ( ; l; l = l->next )
        result = list_append( result, glob_recursive( object_str( l->value ) ) );
    return result;
}


/*
 * builtin_match() - MATCH rule, regexp matching.
 */

LIST * builtin_match( FRAME * frame, int flags )
{
    LIST * l;
    LIST * r;
    LIST * result = 0;

    string buf[ 1 ];
    string_new( buf );

    /* For each pattern */

    for ( l = lol_get( frame->args, 0 ); l; l = l->next )
    {
        /* Result is cached and intentionally never freed. */
        regexp * re = regex_compile( l->value );

        /* For each string to match against. */
        for ( r = lol_get( frame->args, 1 ); r; r = r->next )
        {
            if ( regexec( re, object_str( r->value ) ) )
            {
                int i;
                int top;

                /* Find highest parameter */

                for ( top = NSUBEXP; top-- > 1; )
                    if ( re->startp[ top ] )
                        break;

                /* And add all parameters up to highest onto list. */
                /* Must have parameters to have results! */
                for ( i = 1; i <= top; ++i )
                {
                    string_append_range( buf, re->startp[ i ], re->endp[ i ] );
                    result = list_new( result, object_new( buf->value ) );
                    string_truncate( buf, 0 );
                }
            }
        }
    }

    string_free( buf );
    return result;
}

LIST * builtin_split_by_characters( FRAME * frame, int flags )
{
    LIST * l1 = lol_get( frame->args, 0 );
    LIST * l2 = lol_get( frame->args, 1 );

    LIST * result = L0;
    
    string buf[ 1 ];

    const char * delimiters = object_str( l2->value );
    char * t;

    string_copy( buf, object_str( l1->value ) );

    t = strtok( buf->value, delimiters) ;
    while ( t )
    {
        result = list_new( result, object_new( t ) );
        t = strtok( NULL, delimiters );
    }

    string_free( buf );

    return result;
}

LIST * builtin_hdrmacro( FRAME * frame, int flags )
{
  LIST * l = lol_get( frame->args, 0 );

  for ( ; l; l = list_next( l ) )
  {
    TARGET * t = bindtarget( l->value );

    /* Scan file for header filename macro definitions. */
    if ( DEBUG_HEADER )
        printf( "scanning '%s' for header file macro definitions\n",
            object_str( l->value ) );

    macro_headers( t );
  }

  return L0;
}


/*
 * builtin_rulenames() - RULENAMES ( MODULE ? ).
 *
 * Returns a list of the non-local rule names in the given MODULE. If MODULE is
 * not supplied, returns the list of rule names in the global module.
 */

static void add_rule_name( void * r_, void * result_ )
{
    RULE * r = (RULE *)r_;
    LIST * * result = (LIST * *)result_;
    if ( r->exported )
        *result = list_new( *result, object_copy( r->name ) );
}


LIST * builtin_rulenames( FRAME * frame, int flags )
{
    LIST * arg0 = lol_get( frame->args, 0 );
    LIST * result = L0;
    module_t * source_module = bindmodule( arg0 ? arg0->value : 0 );

    if ( source_module->rules )
        hashenumerate( source_module->rules, add_rule_name, &result );
    return result;
}


/*
 * builtin_varnames() - VARNAMES ( MODULE ? ).
 *
 * Returns a list of the variable names in the given MODULE. If MODULE is not
 * supplied, returns the list of variable names in the global module.
 */

/* helper function for builtin_varnames(), below. Used with hashenumerate, will
 * prepend the key of each element to the list
 */
static void add_hash_key( void * np, void * result_ )
{
    LIST * * result = (LIST * *)result_;
    *result = list_new( *result, object_copy( *(OBJECT * *)np ) );
}


static struct hash * get_running_module_vars()
{
    struct hash * dummy;
    struct hash * vars = NULL;
    /* Get the global variables pointer (that of the currently running module).
     */
    var_hash_swap( &vars );
    dummy = vars;
    /* Put the global variables pointer in its right place. */
    var_hash_swap( &dummy );
    return vars;
}


LIST * builtin_varnames( FRAME * frame, int flags )
{
    LIST * arg0 = lol_get( frame->args, 0 );
    LIST * result = L0;
    module_t * source_module = bindmodule( arg0 ? arg0->value : 0 );

    /* The running module _always_ has its 'variables' member set to NULL due to
     * the way enter_module() and var_hash_swap() work.
     */
    struct hash * vars = source_module == frame->module
        ? get_running_module_vars()
        : source_module->variables;

    if ( vars )
        hashenumerate( vars, add_hash_key, &result );
    return result;
}


/*
 * builtin_delete_module() - MODULE ?.
 *
 * Clears all rules and variables from the given module.
 */

LIST * builtin_delete_module( FRAME * frame, int flags )
{
    LIST     * arg0 = lol_get( frame->args, 0 );
    LIST     * result = L0;
    module_t * source_module = bindmodule( arg0 ? arg0->value : 0 );
    delete_module( source_module );
    return result;
}


static void unknown_rule( FRAME * frame, const char * key, module_t * module, OBJECT * rule_name )
{
    const char * module_name = module->name ? object_str( module->name ) : "";
    backtrace_line( frame->prev );
    if ( module->name )
    {
        printf( "%s error: rule \"%s\" unknown in module \"%s.\"\n", key, object_str( rule_name ), object_str( module->name ) );
    }
    else
    {
        printf( "%s error: rule \"%s\" unknown in module \"\"\n", key, object_str( rule_name ) );
    }
    backtrace( frame->prev );
    exit( 1 );
}


/*
 * builtin_import() - IMPORT
 * (
 *     SOURCE_MODULE ? :
 *     SOURCE_RULES  * :
 *     TARGET_MODULE ? :
 *     TARGET_RULES  * :
 *     LOCALIZE      ?
 * )
 *
 * The IMPORT rule imports rules from the SOURCE_MODULE into the TARGET_MODULE
 * as local rules. If either SOURCE_MODULE or TARGET_MODULE is not supplied, it
 * refers to the global module. SOURCE_RULES specifies which rules from the
 * SOURCE_MODULE to import; TARGET_RULES specifies the names to give those rules
 * in TARGET_MODULE. If SOURCE_RULES contains a name which doesn't correspond to
 * a rule in SOURCE_MODULE, or if it contains a different number of items than
 * TARGET_RULES, an error is issued. If LOCALIZE is specified, the rules will be
 * executed in TARGET_MODULE, with corresponding access to its module local
 * variables.
 */

LIST * builtin_import( FRAME * frame, int flags )
{
    LIST * source_module_list = lol_get( frame->args, 0 );
    LIST * source_rules       = lol_get( frame->args, 1 );
    LIST * target_module_list = lol_get( frame->args, 2 );
    LIST * target_rules       = lol_get( frame->args, 3 );
    LIST * localize           = lol_get( frame->args, 4 );

    module_t * target_module =
        bindmodule( target_module_list ? target_module_list->value : 0 );
    module_t * source_module =
        bindmodule( source_module_list ? source_module_list->value : 0 );

    LIST * source_name;
    LIST * target_name;

    for ( source_name = source_rules, target_name = target_rules;
          source_name && target_name;
          source_name = list_next( source_name ),
          target_name = list_next( target_name ) )
    {
        RULE   r_;
        RULE * r = &r_;
        RULE * imported;
        r_.name = source_name->value;

        if ( !source_module->rules ||
            !hashcheck( source_module->rules, (HASHDATA * *)&r ) )
            unknown_rule( frame, "IMPORT", source_module, r_.name );

        imported = import_rule( r, target_module, target_name->value );
        if ( localize )
            imported->module = target_module;
        /* This rule is really part of some other module. Just refer to it here,
         * but do not let it out.
         */
        imported->exported = 0;
    }

    if ( source_name || target_name )
    {
        backtrace_line( frame->prev );
        printf( "import error: length of source and target rule name lists don't match!\n" );
        printf( "    source: " );
        list_print( source_rules );
        printf( "\n    target: " );
        list_print( target_rules );
        printf( "\n" );
        backtrace( frame->prev );
        exit( 1 );
    }

    return L0;
}


/*
 * builtin_export() - EXPORT ( MODULE ? : RULES * ).
 *
 * The EXPORT rule marks RULES from the SOURCE_MODULE as non-local (and thus
 * exportable). If an element of RULES does not name a rule in MODULE, an error
 * is issued.
 */

LIST * builtin_export( FRAME * frame, int flags )
{
    LIST     * module_list = lol_get( frame->args, 0 );
    LIST     * rules       = lol_get( frame->args, 1 );
    module_t * m           = bindmodule( module_list ? module_list->value : 0 );

    for ( ; rules; rules = list_next( rules ) )
    {
        RULE   r_;
        RULE * r = &r_;
        r_.name = rules->value;

        if ( !m->rules || !hashcheck( m->rules, (HASHDATA * *)&r ) )
            unknown_rule( frame, "EXPORT", m, r_.name );

        r->exported = 1;
    }
    return L0;
}


/*
 * get_source_line() - Retrieve the file and line number that should be
 * indicated for a given procedure in debug output or an error backtrace.
 */

static void get_source_line( FRAME * frame, const char * * file, int * line )
{
    if ( frame->file )
    {
        const char * f = object_str( frame->file );
        int    l = frame->line;
        if ( !strcmp( f, "+" ) )
        {
            f = "jambase.c";
            l += 3;
        }
        *file = f;
        *line = l;
    }
    else
    {
        *file = "(builtin)";
        *line = -1;
    }
}


void print_source_line( FRAME * frame )
{
    const char * file;
    int    line;

    get_source_line( frame, &file, &line );
    if ( line < 0 )
        printf( "(builtin):" );
    else
        printf( "%s:%d:", file, line );
}


/*
 * backtrace_line() - print a single line of error backtrace for the given
 * frame.
 */

void backtrace_line( FRAME * frame )
{
    if ( frame == 0 )
    {
        printf( "(no frame):" );
    }
    else
    {
        print_source_line( frame );
        printf( " in %s\n", frame->rulename );
    }
}


/*
 * backtrace() - Print the entire backtrace from the given frame to the Jambase
 * which invoked it.
 */

void backtrace( FRAME * frame )
{
    if ( !frame ) return;
    while ( ( frame = frame->prev ) )
        backtrace_line( frame );
}


/*
 * builtin_backtrace() - A Jam version of the backtrace function, taking no
 * arguments and returning a list of quadruples: FILENAME LINE MODULE. RULENAME
 * describing each frame. Note that the module-name is always followed by a
 * period.
 */

LIST * builtin_backtrace( FRAME * frame, int flags )
{
    LIST * levels_arg = lol_get( frame->args, 0 );
    int levels = levels_arg ? atoi( object_str( levels_arg->value ) ) : (int)( (unsigned int)(-1) >> 1 ) ;

    LIST * result = L0;
    for ( ; ( frame = frame->prev ) && levels ; --levels )
    {
        const char * file;
        int    line;
        char   buf[32];
        string module_name[1];
        get_source_line( frame, &file, &line );
        sprintf( buf, "%d", line );
        string_new( module_name );
        if ( frame->module->name )
        {
            string_append( module_name, object_str( frame->module->name ) );
            string_append( module_name, "." );
        }
        result = list_new( result, object_new( file ) );
        result = list_new( result, object_new( buf ) );
        result = list_new( result, object_new( module_name->value ) );
        result = list_new( result, object_new( frame->rulename ) );
        string_free( module_name );
    }
    return result;
}


/*
 * builtin_caller_module() - CALLER_MODULE ( levels ? )
 *
 * If levels is not supplied, returns the name of the module of the rule which
 * called the one calling this one. If levels is supplied, it is interpreted as
 * an integer specifying a number of additional levels of call stack to traverse
 * in order to locate the module in question. If no such module exists, returns
 * the empty list. Also returns the empty list when the module in question is
 * the global module. This rule is needed for implementing module import
 * behavior.
 */

LIST * builtin_caller_module( FRAME * frame, int flags )
{
    LIST * levels_arg = lol_get( frame->args, 0 );
    int levels = levels_arg ? atoi( object_str( levels_arg->value ) ) : 0 ;

    int i;
    for ( i = 0; ( i < levels + 2 ) && frame->prev; ++i )
        frame = frame->prev;

    if ( frame->module == root_module() )
        return L0;
    else
        return list_new( L0, object_copy( frame->module->name ) );
}


/*
 * Return the current working directory.
 *
 * Usage: pwd = [ PWD ] ;
 */

LIST * builtin_pwd( FRAME * frame, int flags )
{
    return pwd();
}


/*
 * Adds targets to the list of target that jam will attempt to update.
 */

LIST * builtin_update( FRAME * frame, int flags )
{
    LIST * result = list_copy( L0, targets_to_update() );
    LIST * arg1 = lol_get( frame->args, 0 );
    clear_targets_to_update();
    for ( ; arg1; arg1 = list_next( arg1 ) )
        mark_target_for_updating( object_copy( arg1->value ) );
    return result;
}

extern int anyhow;
int last_update_now_status;

/* Takes a list of target names as first argument, and immediately
   updates them.
   Second parameter, if specified, if the descriptor (converted to a string)
   of a log file where all build output is redirected.
   Third parameter, if non-empty, specifies that the -n option should have
   no effect -- that is, all out-of-date targets should be rebuild.
*/
LIST * builtin_update_now( FRAME * frame, int flags )
{
    LIST * targets = lol_get( frame->args, 0 );
    LIST * log = lol_get( frame->args, 1 );
    LIST * force = lol_get( frame->args, 2 );
    LIST * continue_ = lol_get( frame->args, 3 );
    int status = 0;
    int original_stdout;
    int original_stderr;
    int n;
    int targets_count;
    OBJECT * * targets2;
    int i;
    int original_noexec;
    int original_quitquick;
	

    if ( log )
    {
        int fd = atoi( object_str( log->value ) );
        /* Redirect stdout and stderr, temporary, to the log file.  */
        original_stdout = dup( 0 );
        original_stderr = dup( 1 );
        dup2 ( fd, 0 );
        dup2 ( fd, 1 );
    }

    if ( force )
    {
        original_noexec = globs.noexec;
        globs.noexec = 0;
        original_quitquick = globs.quitquick;
        globs.quitquick = 0;
    }

    if ( continue_ )
    {
        original_quitquick = globs.quitquick;
        globs.quitquick = 0;
    }

    targets_count = list_length( targets );
    targets2 = (OBJECT * *)BJAM_MALLOC( targets_count * sizeof( OBJECT * ) );    
    for (i = 0 ; targets; targets = list_next( targets ) )
        targets2[ i++ ] = targets->value;
    status |= make( targets_count, targets2, anyhow);
    BJAM_FREE( (void *)targets2 );

    if (force)
    {
        globs.noexec = original_noexec;
        globs.quitquick = original_quitquick;
    }

    if ( continue_ )
    {
        globs.quitquick = original_quitquick;
    }

    if ( log )
    {
        /* Flush whatever stdio might have buffered, while descriptions
           0 and 1 still refer to the log file.  */
        fflush( stdout );
        fflush( stderr );
        dup2( original_stdout, 0 );
        dup2( original_stderr, 1 );
        close( original_stdout );
        close( original_stderr );
    }

    last_update_now_status = status;
	
    if ( status == 0 )
        return list_new( L0, object_new( "ok" ) );
    else
        return L0;
}

LIST * builtin_search_for_target( FRAME * frame, int flags )
{
    LIST * arg1 = lol_get( frame->args, 0 );
    LIST * arg2 = lol_get( frame->args, 1 );
    TARGET * t = search_for_target( arg1->value, arg2 );
    return list_new( L0, object_copy( t->name ) );
}


LIST * builtin_import_module( FRAME * frame, int flags )
{
    LIST * arg1 = lol_get( frame->args, 0 );
    LIST * arg2 = lol_get( frame->args, 1 );
    module_t * m = arg2 ? bindmodule( arg2->value ) : root_module();
    import_module( arg1, m );
    return L0;
}


LIST * builtin_imported_modules( FRAME * frame, int flags )
{
    LIST * arg0 = lol_get( frame->args, 0 );
    return imported_modules( bindmodule( arg0 ? arg0->value : 0 ) );
}


LIST * builtin_instance( FRAME * frame, int flags )
{
    LIST * arg1 = lol_get( frame->args, 0 );
    LIST * arg2 = lol_get( frame->args, 1 );
    module_t * const instance     = bindmodule( arg1->value );
    module_t * const class_module = bindmodule( arg2->value );
    instance->class_module = class_module;
    return L0;
}


LIST * builtin_sort( FRAME * frame, int flags )
{
    LIST * arg1 = lol_get( frame->args, 0 );
    return list_sort( arg1 );
}


LIST * builtin_normalize_path( FRAME * frame, int flags )
{
    LIST * arg = lol_get( frame->args, 0 );

    /* First, we iterate over all '/'-separated elements, starting from the end
     * of string. If we see a '..', we remove a previous path elements. If we
     * see '.', we remove it. The removal is done by overwriting data using '\1'
     * in the string. After the whole string has been processed, we do a second
     * pass, removing all the entered '\1' characters.
     */

    string   in[ 1 ];
    string   out[ 1 ];
	/* Last character of the part of string still to be processed. */
    char   * end;
	/* Working pointer. */
    char   * current;
	/* Number of '..' elements seen and not processed yet. */
    int      dotdots = 0;
    int      rooted  = 0;
    OBJECT * result  = 0;

    /* Make a copy of input: we should not change it. Prepend a '/' before it as
     * a guard for the algorithm later on and remember whether it was originally
     * rooted or not.
     */
    string_new( in );
    string_push_back( in, '/' );
    for ( ; arg; arg = list_next( arg ) )
    {
        if ( object_str( arg->value )[ 0 ] != '\0' )
        {
            if ( in->size == 1 )
                rooted = ( ( object_str( arg->value )[ 0 ] == '/'  ) ||
                           ( object_str( arg->value )[ 0 ] == '\\' ) );
            else
                string_append( in, "/" );
            string_append( in, object_str( arg->value ) );
        }
    }

    /* Convert \ into /. On Windows, paths using / and \ are equivalent, and we
     * want this function to obtain a canonic representation.
     */
    for ( current = in->value, end = in->value + in->size;
        current < end; ++current )
        if ( *current == '\\' )
            *current = '/';

    /* Now we remove any extra path elements by overwriting them with '\1'
     * characters and cound how many more unused '..' path elements there are
     * remaining. Note that each remaining path element with always starts with
     * a '/' character.
     */
    for ( end = in->value + in->size - 1; end >= in->value; )
    {
        /* Set 'current' to the next occurence of '/', which always exists. */
        for ( current = end; *current != '/'; --current );

        if ( current == end )
        {
            /* Found a trailing or duplicate '/'. Remove it. */
            *current = '\1';
        }
        else if ( ( end - current == 1 ) && ( *(current + 1) == '.' ) )
        {
            /* Found '/.'. Remove them all. */
            *current = '\1';
            *(current + 1) = '\1';
        }
        else if ( ( end - current == 2 ) && ( *(current + 1) == '.' ) && ( *(current + 2) == '.' ) )
        {
            /* Found '/..'. Remove them all. */
            *current = '\1';
            *(current + 1) = '\1';
            *(current + 2) = '\1';
            ++dotdots;
        }
        else if ( dotdots )
        {
            memset( current, '\1', end - current + 1 );
            --dotdots;
        }
        end = current - 1;
    }

    string_new( out );

    /* Now we know that we need to add exactly dotdots '..' path elements to the
     * front and that our string is either empty or has a '/' as its first
     * significant character. If we have any dotdots remaining then the passed
     * path must not have been rooted or else it is invalid we return an empty
     * list.
     */
    if ( dotdots )
    {
        if ( rooted ) return L0;
        do
            string_append( out, "/.." );
        while ( --dotdots );
    }

    /* Now we actually remove all the path characters marked for removal. */
    for ( current = in->value; *current; ++current )
        if ( *current != '\1' )
            string_push_back( out, *current );

    /* Here we know that our string contains no '\1' characters and is either
     * empty or has a '/' as its initial character. If the original path was not
     * rooted and we have a non-empty path we need to drop the initial '/'. If
     * the original path was rooted and we have an empty path we need to add
     * back the '/'.
     */
    result = object_new( out->size ? out->value + !rooted : ( rooted ? "/" : "." ) );

    string_free( out );
    string_free( in );

    return list_new( 0, result );
}


LIST * builtin_native_rule( FRAME * frame, int flags )
{
    LIST * module_name = lol_get( frame->args, 0 );
    LIST * rule_name = lol_get( frame->args, 1 );

    module_t * module = bindmodule( module_name->value );

    native_rule_t n;
    native_rule_t * np = &n;
    n.name = rule_name->value;
    if ( module->native_rules && hashcheck( module->native_rules, (HASHDATA * *)&np ) )
    {
        args_refer( np->arguments );
        new_rule_body( module, np->name, np->arguments, np->procedure, 1 );
    }
    else
    {
        backtrace_line( frame->prev );
        printf( "error: no native rule \"%s\" defined in module \"%s.\"\n",
                object_str( n.name ), object_str( module->name ) );
        backtrace( frame->prev );
        exit( 1 );
    }
    return L0;
}


LIST * builtin_has_native_rule( FRAME * frame, int flags )
{
    LIST * module_name = lol_get( frame->args, 0 );
    LIST * rule_name   = lol_get( frame->args, 1 );
    LIST * version     = lol_get( frame->args, 2 );

    module_t * module = bindmodule( module_name->value );

    native_rule_t n;
    native_rule_t * np = &n;
    n.name = rule_name->value;
    if ( module->native_rules && hashcheck( module->native_rules, (HASHDATA * *)&np ) )
    {
        int expected_version = atoi( object_str( version->value ) );
        if ( np->version == expected_version )
            return list_new( 0, object_new( "true" ) );
    }
    return L0;
}


LIST * builtin_user_module( FRAME * frame, int flags )
{
    LIST * module_name = lol_get( frame->args, 0 );
    for ( ; module_name; module_name = module_name->next )
    {
        module_t * m = bindmodule( module_name->value );
        m->user_module = 1;
    }
    return L0;
}


LIST * builtin_nearest_user_location( FRAME * frame, int flags )
{
    FRAME * nearest_user_frame =
        frame->module->user_module ? frame : frame->prev_user;
    if ( !nearest_user_frame )
        return L0;

    {
        LIST * result = 0;
        const char * file;
        int    line;
        char   buf[32];

        get_source_line( nearest_user_frame, &file, &line );
        sprintf( buf, "%d", line );
        result = list_new( result, object_new( file ) );
        result = list_new( result, object_new( buf ) );
        return result;
    }
}


LIST * builtin_check_if_file( FRAME * frame, int flags )
{
    LIST * name = lol_get( frame->args, 0 );
    return file_is_file( name->value ) == 1
        ? list_new( 0, object_new( "true" ) )
        : L0 ;
}


LIST * builtin_md5( FRAME * frame, int flags )
{
    LIST * l = lol_get( frame->args, 0 );
    const char* s = object_str( l->value );

    md5_state_t state;
    md5_byte_t digest[16];
    char hex_output[16*2 + 1];

    int di;

    md5_init( &state );
    md5_append( &state, (const md5_byte_t *)s, strlen(s) );
    md5_finish( &state, digest );

    for (di = 0; di < 16; ++di)
        sprintf( hex_output + di * 2, "%02x", digest[di] );

    return list_new( L0, object_new( hex_output ) );
}

LIST *builtin_file_open( FRAME * frame, int flags )
{
    const char * name = object_str( lol_get( frame->args, 0 )->value );
    const char * mode = object_str( lol_get( frame->args, 1 )->value );
    int fd;
    char buffer[sizeof("4294967295")];

    if ( strcmp(mode, "w") == 0 )
    {
        fd = open( name, O_WRONLY|O_CREAT|O_TRUNC, 0666 );
    }
    else
    {
        fd = open( name, O_RDONLY );
    }

    if (fd != -1)
    {
        sprintf( buffer, "%d", fd );
        return list_new( L0, object_new( buffer ) );
    }
    else
    {
        return L0;
    }
}

LIST *builtin_pad( FRAME * frame, int flags )
{
    OBJECT * string = lol_get( frame->args, 0 )->value;
    const char * width_s = object_str( lol_get( frame->args, 1 )->value );

    int current = strlen( object_str( string ) );
    int desired = atoi( width_s );
    if (current >= desired)
        return list_new (L0, object_copy( string ) );
    else
    {
        char * buffer = BJAM_MALLOC( desired + 1 );
        int i;
        LIST * result;

        strcpy( buffer, object_str( string ) );
        for ( i = current; i < desired; ++i )
            buffer[i] = ' ';
        buffer[desired] = '\0';
        result = list_new( L0, object_new( buffer ) );
        BJAM_FREE( buffer );
        return result;
    }
}

LIST *builtin_precious( FRAME * frame, int flags )
{
    LIST * targets = lol_get(frame->args, 0);

    for ( ; targets; targets = list_next( targets ) )    
    {
        TARGET* t = bindtarget( targets->value );
        t->flags |= T_FLAG_PRECIOUS;
    }

    return L0;
}

LIST *builtin_self_path( FRAME * frame, int flags )
{
    extern const char * saved_argv0;
    char * p = executable_path( saved_argv0 );
    if ( p )
    {
        LIST* result = list_new( 0, object_new( p ) );
        free( p );
        return result;
    }
    else
    {
        return L0;
    }
}

LIST *builtin_makedir( FRAME * frame, int flags )
{
    LIST * path = lol_get( frame->args, 0 );

    if ( file_mkdir( object_str( path->value ) ) == 0 )
    {
        LIST * result = list_new ( L0, object_copy( path->value ) );
        return result;
    }
    else
    {
        return L0;
    }    
}

LIST * builtin_rescan( FRAME * frame, int flags )
{
    LIST * directory = lol_get( frame->args, 0 );

    if ( !directory )
    {
        /* clear file and timestamp hashes */
        file_free_all();
        time_free_all();
    }
    else
    {
        for ( ; directory; directory = list_next( directory ) )
        {
            file_free( directory->value, 0 );
            time_free( directory->value );
        }
    }
    return L0;
}

#ifdef HAVE_PYTHON

LIST * builtin_python_import_rule( FRAME * frame, int flags )
{
    static int first_time = 1;
    const char * python_module   = object_str( lol_get( frame->args, 0 )->value );
    const char * python_function = object_str( lol_get( frame->args, 1 )->value );
    OBJECT     * jam_module      = lol_get( frame->args, 2 )->value;
    OBJECT     * jam_rule        = lol_get( frame->args, 3 )->value;

    PyObject * pName;
    PyObject * pModule;
    PyObject * pDict;
    PyObject * pFunc;

    if ( first_time )
    {
        /* At the first invocation, we add the value of the global
         * EXTRA_PYTHONPATH to the sys.path Python variable.
         */
        LIST * extra = 0;
        module_t * outer_module = frame->module;

        first_time = 0;

        if ( outer_module != root_module() )
        {
            exit_module( outer_module );
            enter_module( root_module() );
        }

        extra = var_get( constant_extra_pythonpath );

        if ( outer_module != root_module() )
        {
             exit_module( root_module() );
             enter_module( outer_module );
        }

        for ( ; extra; extra = extra->next )
        {
            string buf[ 1 ];
            string_new( buf );
            string_append( buf, "import sys\nsys.path.append(\"" );
            string_append( buf, object_str( extra->value ) );
            string_append( buf, "\")\n" );
            PyRun_SimpleString( buf->value );
            string_free( buf );
        }
    }

    pName   = PyString_FromString( python_module );
    pModule = PyImport_Import( pName );
    Py_DECREF( pName );

    if ( pModule != NULL )
    {
        pDict = PyModule_GetDict( pModule );
        pFunc = PyDict_GetItemString( pDict, python_function );

        if ( pFunc && PyCallable_Check( pFunc ) )
        {
            module_t * m = bindmodule( jam_module );
            RULE * r = bindrule( jam_rule, m );

            /* Make pFunc owned. */
            Py_INCREF( pFunc );

            r->python_function = pFunc;
        }
        else
        {
            if ( PyErr_Occurred() )
                PyErr_Print();
            fprintf( stderr, "Cannot find function \"%s\"\n", python_function );
        }
        Py_DECREF( pModule );
    }
    else
    {
        PyErr_Print();
        fprintf( stderr, "Failed to load \"%s\"\n", python_module );
    }
    return L0;

}

#endif

void lol_build( LOL * lol, const char * * elements )
{
    LIST * l = L0;
    lol_init( lol );

    while ( elements && *elements )
    {
        if ( !strcmp( *elements, ":" ) )
        {
            lol_add( lol, l );
            l = L0 ;
        }
        else
        {
            l = list_new( l, object_new( *elements ) );
        }
        ++elements;
    }

    if ( l != L0 )
        lol_add( lol, l );
}


#ifdef HAVE_PYTHON

/*
 * Calls the bjam rule specified by name passed in 'args'. The name is looked up
 * in the context of bjam's 'python_interface' module. Returns the list of
 * string retured by the rule.
 */

PyObject* bjam_call( PyObject * self, PyObject * args )
{
    FRAME    inner[ 1 ];
    LIST   * result;
    PARSE  * p;
    OBJECT * rulename;

    /* Build up the list of arg lists. */
    frame_init( inner );
    inner->prev = 0;
    inner->prev_user = 0;
    inner->module = bindmodule( constant_python_interface );

    /* Extract the rule name and arguments from 'args'. */

    /* PyTuple_GetItem returns borrowed reference. */
    rulename = object_new( PyString_AsString( PyTuple_GetItem( args, 0 ) ) );
    {
        int i = 1;
        int size = PyTuple_Size( args );
        for ( ; i < size; ++i )
        {
            PyObject * a = PyTuple_GetItem( args, i );
            if ( PyString_Check( a ) )
            {
                lol_add( inner->args, list_new( 0, object_new(
                    PyString_AsString( a ) ) ) );
            }
            else if ( PySequence_Check( a ) )
            {
                LIST * l = 0;
                int s = PySequence_Size( a );
                int i = 0;
                for ( ; i < s; ++i )
                {
                    /* PySequence_GetItem returns new reference. */
                    PyObject * e = PySequence_GetItem( a, i );
                    char * s = PyString_AsString( e );
                    if ( !s )
                    {
                        printf( "Invalid parameter type passed from Python\n" );
                        exit( 1 );
                    }
                    l = list_new( l, object_new( s ) );
                    Py_DECREF( e );
                }
                lol_add( inner->args, l );
            }
        }
    }

    result = evaluate_rule( rulename, inner );
    object_free( rulename );

    frame_free( inner );

    /* Convert the bjam list into a Python list result. */
    {
        PyObject * pyResult = PyList_New( list_length( result ) );
        int i = 0;
        while ( result )
        {
            PyList_SetItem( pyResult, i, PyString_FromString( object_str( result->value ) ) );
            result = list_next( result );
            i += 1;
        }
        list_free( result );
        return pyResult;
    }
}


/*
 * Accepts four arguments: 
 * - module name
 * - rule name,
 * - Python callable. 
 * - (optional) bjam language function signature.
 * Creates a bjam rule with the specified name in the specified module, which will
 * invoke the Python callable.
 */

PyObject * bjam_import_rule( PyObject * self, PyObject * args )
{
    char     * module;
    char     * rule;
    PyObject * func;
    PyObject * bjam_signature = NULL;
    module_t * m;
    RULE     * r;
    OBJECT   * module_name;
    OBJECT   * rule_name;

    if ( !PyArg_ParseTuple( args, "ssO|O:import_rule", 
                            &module, &rule, &func, &bjam_signature ) )
        return NULL;

    if ( !PyCallable_Check( func ) )
    {
        PyErr_SetString( PyExc_RuntimeError,
                        "Non-callable object passed to bjam.import_rule" );
        return NULL;
    }

    module_name = *module ? object_new( module ) : 0;
    m = bindmodule( module_name );
    if( module_name )
    {
        object_free( module_name );
    }
    rule_name = object_new( rule );
    r = bindrule( rule_name, m );
    object_free( rule_name );

    /* Make pFunc owned. */
    Py_INCREF( func );

    r->python_function = func;
    r->arguments = 0;

    if (bjam_signature)
    {
        argument_list * arg_list = args_new();
        Py_ssize_t i;

        Py_ssize_t s = PySequence_Size (bjam_signature);
        for (i = 0; i < s; ++i)
        {
            PyObject* v = PySequence_GetItem (bjam_signature, i);
            lol_add(arg_list->data, list_from_python (v));
            Py_DECREF(v);
        }
        r->arguments = arg_list;
    }

    Py_INCREF( Py_None );
    return Py_None;
}


/*
 * Accepts four arguments:
 *  - an action name
 *  - an action body
 *  - a list of variable that will be bound inside the action
 *  - integer flags.
 *  Defines an action on bjam side.
 */

PyObject * bjam_define_action( PyObject * self, PyObject * args )
{
    char     * name;
    char     * body;
    module_t * m;
    PyObject * bindlist_python;
    int        flags;
    LIST     * bindlist = L0;
    int        n;
    int        i;
    OBJECT   * name_str;
    OBJECT   * body_str;

    if ( !PyArg_ParseTuple( args, "ssO!i:define_action", &name, &body,
                          &PyList_Type, &bindlist_python, &flags ) )
        return NULL;

    n = PyList_Size( bindlist_python );
    for ( i = 0; i < n; ++i )
    {
        PyObject * next = PyList_GetItem( bindlist_python, i );
        if ( !PyString_Check( next ) )
        {
            PyErr_SetString( PyExc_RuntimeError,
                            "bind list has non-string type" );
            return NULL;
        }
        bindlist = list_new( bindlist, object_new( PyString_AsString( next ) ) );
    }

    name_str = object_new( name );
    body_str = object_new( body );
    new_rule_actions( root_module(), name_str, body_str, bindlist, flags );
    object_free( body_str );
    object_free( name_str );

    Py_INCREF( Py_None );
    return Py_None;
}


/*
 * Returns the value of a variable in root Jam module.
 */

PyObject * bjam_variable( PyObject * self, PyObject * args )
{
    char     * name;
    LIST     * value;
    PyObject * result;
    int        i;
    OBJECT   * varname;

    if ( !PyArg_ParseTuple( args, "s", &name ) )
        return NULL;

    enter_module( root_module() );
    varname = object_new( name );
    value = var_get( varname );
    object_free( varname );
    exit_module( root_module() );

    result = PyList_New( list_length( value ) );
    for ( i = 0; value; value = list_next( value ), ++i )
        PyList_SetItem( result, i, PyString_FromString( object_str( value->value ) ) );

    return result;
}


PyObject * bjam_backtrace( PyObject * self, PyObject * args )
{
    PyObject     * result = PyList_New( 0 );
    struct frame * f = frame_before_python_call;

    for ( ; f = f->prev; )
    {
        PyObject   * tuple = PyTuple_New( 4 );
        const char * file;
        int          line;
        char         buf[ 32 ];
        string module_name[1];

        get_source_line( f, &file, &line );
        sprintf( buf, "%d", line );
        string_new( module_name );
        if ( f->module->name )
        {
            string_append( module_name, object_str( f->module->name ) );
            string_append( module_name, "." );
        }

        /* PyTuple_SetItem steals reference. */
        PyTuple_SetItem( tuple, 0, PyString_FromString( file            ) );
        PyTuple_SetItem( tuple, 1, PyString_FromString( buf             ) );
        PyTuple_SetItem( tuple, 2, PyString_FromString( module_name->value ) );
        PyTuple_SetItem( tuple, 3, PyString_FromString( f->rulename ) );

        string_free( module_name );

        PyList_Append( result, tuple );
        Py_DECREF( tuple );
    }
    return result;
}

PyObject * bjam_caller( PyObject * self, PyObject * args )
{
    const char * s =  frame_before_python_call->prev->module->name ?
        object_str( frame_before_python_call->prev->module->name ) :
        "";
    return PyString_FromString( s );
}

#endif  /* #ifdef HAVE_PYTHON */


#ifdef HAVE_POPEN

#if defined(_MSC_VER) || defined(__BORLANDC__)
    #define popen windows_popen_wrapper
    #define pclose _pclose

    /*
     * This wrapper is a workaround for a funny _popen() feature on Windows
     * where it eats external quotes in some cases. The bug seems to be related
     * to the quote stripping functionality used by the Windows cmd.exe
     * interpreter when its /S is not specified.
     *
     * Cleaned up quote from the cmd.exe help screen as displayed on Windows XP
     * SP3:
     *
     *   1. If all of the following conditions are met, then quote characters on
     *      the command line are preserved:
     *
     *       - no /S switch
     *       - exactly two quote characters
     *       - no special characters between the two quote characters, where
     *         special is one of: &<>()@^|
     *       - there are one or more whitespace characters between the two quote
     *         characters
     *       - the string between the two quote characters is the name of an
     *         executable file.
     *
     *   2. Otherwise, old behavior is to see if the first character is a quote
     *      character and if so, strip the leading character and remove the last
     *      quote character on the command line, preserving any text after the
     *      last quote character.
     *
     * This causes some commands containing quotes not to be executed correctly.
     * For example:
     *
     *   "\Long folder name\aaa.exe" --name="Jurko" --no-surname
     *
     * would get its outermost quotes stripped and would be executed as:
     *
     *   \Long folder name\aaa.exe" --name="Jurko --no-surname
     *
     * which would report an error about '\Long' not being a valid command.
     *
     * cmd.exe help seems to indicate it would be enough to add an extra space
     * character in front of the command to avoid this but this does not work,
     * most likely due to the shell first stripping all leading whitespace
     * characters from the command.
     *
     * Solution implemented here is to quote the whole command in case it
     * contains any quote characters. Note thought this will not work correctly
     * should Windows ever 'fix' this feature.
     *                                               (03.06.2008.) (Jurko)
     */
    static FILE * windows_popen_wrapper( const char * command, const char * mode )
    {
        int extra_command_quotes_needed = ( strchr( command, '"' ) != 0 );
        string quoted_command;
        FILE * result;

        if ( extra_command_quotes_needed )
        {
            string_new( &quoted_command );
            string_append( &quoted_command, "\"" );
            string_append( &quoted_command, command );
            string_append( &quoted_command, "\"" );
            command = quoted_command.value;
        }

        result = _popen( command, "r" );

        if ( extra_command_quotes_needed )
            string_free( &quoted_command );

        return result;
    }
#endif


static char * rtrim( char * s )
{
    char * p = s;
    while ( *p ) ++p;
    for ( --p; p >= s && isspace( *p ); *p-- = 0 );
    return s;
}

LIST * builtin_shell( FRAME * frame, int flags )
{
    LIST   * command = lol_get( frame->args, 0 );
    LIST   * result = 0;
    string   s;
    int      ret;
    char     buffer[ 1024 ];
    FILE   * p = NULL;
    int      exit_status = -1;
    int      exit_status_opt = 0;
    int      no_output_opt = 0;
    int      strip_eol_opt = 0;

    /* Process the variable args options. */
    {
        int a = 1;
        LIST * arg = lol_get( frame->args, a );
        while ( arg )
        {
            if ( strcmp( "exit-status", object_str( arg->value ) ) == 0 )
            {
                exit_status_opt = 1;
            }
            else if ( strcmp( "no-output", object_str( arg->value ) ) == 0 )
            {
                no_output_opt = 1;
            }
            else if ( strcmp("strip-eol", object_str( arg->value ) ) == 0 )
            {
                strip_eol_opt = 1;
            }
            arg = lol_get( frame->args, ++a );
        }
    }

    /* The following fflush() call seems to be indicated as a workaround for a
     * popen() bug on POSIX implementations related to synhronizing input
     * stream positions for the called and the calling process.
     */
    fflush( NULL );

    p = popen( object_str( command->value ), "r" );
    if ( p == NULL )
        return L0;

    string_new( &s );

    while ( ( ret = fread( buffer, sizeof( char ), sizeof( buffer ) - 1, p ) ) > 0 )
    {
        buffer[ret] = 0;
        if ( !no_output_opt )
        {
            if ( strip_eol_opt )
                rtrim(buffer);
            string_append( &s, buffer );
        }
    }

    exit_status = pclose( p );

    /* The command output is returned first. */
    result = list_new( L0, object_new( s.value ) );
    string_free( &s );

    /* The command exit result next. */
    if ( exit_status_opt )
    {
        if ( WIFEXITED(exit_status) )
            exit_status = WEXITSTATUS(exit_status);
        else
            exit_status = -1;
        sprintf( buffer, "%d", exit_status );
        result = list_new( result, object_new( buffer ) );
    }

    return result;
}

#else  /* #ifdef HAVE_POPEN */

LIST * builtin_shell( FRAME * frame, int flags )
{
    return L0;
}

#endif /* #ifdef HAVE_POPEN */
