/*
 * Copyright (c) 2011-2013 Michael Pippig
 *
 * This file is part of PFFT.
 *
 * PFFT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PFFT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PFFT.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include <stdlib.h>     /* atoi */
#include <stddef.h>     /* ptrdiff_t */
#include <string.h>

#include <pfft.h>
#include <ipfft.h>

static int read_one_arg(
    char **argv, int argNum, unsigned type,
    void *parameter);

static int count_given_args(
    char **argv);
static int enough_args_given(
    char **argv, int neededArgs);
static int strings_match(
    const char *str1, const char *str2);


void PX(get_args)(
    int argc, char **argv, const char *name, int neededArgs, unsigned type,
    void *parameter
    )
{
  int err=0;
  
  for(int iter=0; iter<argc; iter++){
    if ( strings_match( argv[iter], name ) ){
      if ( !enough_args_given(argv+iter, neededArgs) )
        PX(printf)(MPI_COMM_WORLD, "!!! Warning: Not enough command line arguments for %s !!!\n", name);
      else{
	for(int t=0; t<neededArgs || err; t++, iter++)
	  err = read_one_arg(argv+iter+1, t, type, parameter);
	if(err)
	  PX(fprintf)(MPI_COMM_WORLD, stderr, "!!! Error: PFFT_DATATYPE of %s not supported. !!!\n", name);
      }
      return;
    }
  }
}


static int read_one_arg(
    char **argv, int argNum, unsigned type,
    void *parameter
    )
{
  int rtn = 0;
  
  switch(type){
    case PFFT_INT:
      ((int*)parameter)[argNum] = atoi(argv[0]); break;
    case PFFT_PTRDIFF_T:
      ((ptrdiff_t*)parameter)[argNum] = atoi(argv[0]); break;
    case PFFT_FLOAT:
      ((float*)parameter)[argNum] = (float) atof(argv[0]); break;
    case PFFT_DOUBLE:
      ((double*)parameter)[argNum] = atof(argv[0]); break;
    case PFFT_UNSIGNED:
      ((unsigned*)parameter)[argNum] = atoi(argv[0]); break;
    default:
      rtn = 1;
  }
    
    return rtn;
}

static int count_given_args(
    char **argv
    )
{
  int givenArgs;
  for(givenArgs=0; argv[givenArgs+1] != NULL; givenArgs++)
    if (argv[givenArgs+1][0] == '-')
      break;

  return givenArgs;
}

static int enough_args_given(
    char **current_argv, int neededArgs
    )
{
  int givenArgs;
  givenArgs = count_given_args(current_argv);
  return givenArgs >= neededArgs;
}

static int strings_match(
    const char *str1, const char *str2
    )
{
  return !strcmp(str1, str2);
}


