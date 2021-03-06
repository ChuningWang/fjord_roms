/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2017 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Estuary with Sediment Transport Test.
**
** Application flag:   ESTUARY_TEST
** Input script:       ocean_estuary_test.in
**                     sediment_estuary_test.in
*/

/* general */
#define MASKING
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS

/* outputs */
#undef NO_WRITE_GRID
#undef NO_HIS

/* iceplume */
#define ICEPLUME

/* advection, dissipation, pressure grad, etc. */
#define UV_ADV
#define UV_COR

#define TS_U3HADVECTION
#ifdef SOLVE3D
# define TS_C4VADVECTION
#endif

#define UV_VIS2
#define MIX_S_UV
#define TS_DIF2
#undef MIX_S_TS
#define MIX_GEO_TS

#ifdef SOLVE3D
# define DJ_GRADPS
#endif

/* vertical mixing */
#ifdef SOLVE3D
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define RI_SPLINES
#endif

#ifdef SOLVE3D
# define GLS_MIXING
# undef MY25_MIXING

# if defined GLS_MIXING || defined MY25_MIXING
#  define N2S2_HORAVG
#  define CRAIG_BANNER
#  define KANTHA_CLAYSON
#  define CHARNOK
# endif
#endif

#define ANA_DRAG
#define UV_DRAG_GRID
#define UV_LDRAG
#undef UV_QDRAG

#define VISC_GRID
#define DIFF_GRID

/* tides */
#define LTIDES
#ifdef LTIDES
# define SSH_TIDES
# define UV_TIDES
# define ADD_FSOBC
# define ADD_M2OBC
# define RAMP_TIDES
#endif

/* tracers */
#define T_PASSIVE
#define ANA_PASSIVE

/* sediment */
#ifdef SEDIMENT
# define SUSPLOAD
#endif

/* analytical functionals */
#undef ANA_GRID
#define ANA_INITIAL
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_M3OBC
#define ANA_TOBC
#undef ANA_PSOURCE
#ifdef SEDIMENT
# define ANA_SEDIMENT
#endif
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_SPFLUX
#define ANA_BPFLUX
