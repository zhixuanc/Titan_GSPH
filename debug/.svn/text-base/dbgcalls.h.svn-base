#ifndef DBGCALLS__H
#define DBGCALLS__H

typedef struct {
  char filename[15];
  bool write2file;
  bool trace_enabled;
} DbgVars;

void init_debug_vars(int myid, DbgVars *dbstruct);
void dbg_printf(const char *, ...);

#endif // DBGCALLS__H
