/* the following stub is due to the bug in MS_MPI - https://github.com/microsoft/Microsoft-MPI/issues/7
 * this workaround was proposed by Maik Riechert
 */
void __guard_check_icall_fptr(unsigned long ptr) { }