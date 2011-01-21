/* #undef PECOS_BUILD_SHARED_LIBS */
#if defined (_WIN32) && defined (PECOS_BUILD_SHARED_LIBS)
#  if defined(PECOS_LIB_EXPORTS_MODE)
#    define PECOS_LIB_DLL_EXPORT __declspec(dllexport)
#  else
#    define PECOS_LIB_DLL_EXPORT __declspec(dllimport)
#  endif
#else
#  define PECOS_LIB_DLL_EXPORT
#endif
