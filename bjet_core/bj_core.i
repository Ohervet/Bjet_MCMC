%module(packages="bjet_core.bj_core") bj_core
%{
#define SWIG_FILE_WITH_INIT
#include "bj_core.h"
#include "processes_supp_core.h"
%}
%include bj_core.h
%include processes_supp_core.h