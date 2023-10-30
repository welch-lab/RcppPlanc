#include "hw_detect.h"
#include "hwloc.h"

unsigned int get_l1_data_cache(void) {
    hwloc_topology_t topo;
    hwloc_topology_init(&topo);
    hwloc_topology_load(topo);
    hwloc_obj_t myL1 = hwloc_get_obj_by_type(topo, HWLOC_OBJ_L1CACHE, 0);
    hwloc_uint64_t myL1size = myL1->attr->cache.size;
    hwloc_topology_destroy(topo);
    return myL1size;
}

unsigned int get_l2_data_cache(void)
{
    hwloc_topology_t topo;
    hwloc_topology_init(&topo);
    hwloc_topology_load(topo);
    hwloc_obj_t myL2 = hwloc_get_obj_by_type(topo, HWLOC_OBJ_L2CACHE, 0);
    hwloc_uint64_t myL2size = myL2->attr->cache.size;
    hwloc_topology_destroy(topo);
    return myL2size;
}
unsigned int get_num_bound_threads(void) {
    hwloc_topology_t topo;
    hwloc_topology_init(&topo);
    hwloc_topology_load(topo);
    hwloc_cpuset_t cpuset = hwloc_bitmap_alloc();
    hwloc_get_cpubind(topo, cpuset, HWLOC_CPUBIND_PROCESS);
    unsigned int cpuid = 0, numPU = 0;
    hwloc_bitmap_foreach_begin(cpuid, cpuset);
    ++numPU;
    hwloc_bitmap_foreach_end();
    hwloc_bitmap_free(cpuset);
    hwloc_topology_destroy(topo);
    return numPU;
}
