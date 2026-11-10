/* LD_PRELOAD shim: make sync() a no-op. On this WSL2 machine a host-backed filesystem sometimes
   stops answering flush requests, and every sync() then blocks forever (SPAdes calls sync()).
   Assemblers don't need sync() for correctness; the kernel writes the data back regardless. */
void sync(void) {}
