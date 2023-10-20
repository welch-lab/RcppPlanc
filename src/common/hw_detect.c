#include "hw_detect.h"
#include "hwloc.h"

unsigned int get_l1_data_cache() {

    hwloc_topology_t* topo = new hwloc_topology_t;
    hwloc_topology_init(topo);
    hwloc_topology_load(*topo);
    hwloc_obj_t myL1 = hwloc_get_obj_by_type(*topo, HWLOC_OBJ_L1CACHE, 0);
    free(topo);
    return myL1->attr->cache.size;
}

unsigned int get_l2_data_cache()
{

    hwloc_topology_t *topo = new hwloc_topology_t;
    hwloc_topology_init(topo);
    hwloc_topology_load(*topo);
    hwloc_obj_t myL2 = hwloc_get_obj_by_type(*topo, HWLOC_OBJ_L2CACHE, 0);
    free(topo);
    return myL2->attr->cache.size;
}
