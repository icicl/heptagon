#include <stdint.h>


#ifdef __cplusplus
extern "C" {
#endif
    typedef struct Node Node;
    struct Node {
        uint32_t id;
        uint32_t left;
        uint32_t right;
        uint32_t parent;
        double dist;
    };

    void print_node(Node* node_data, uint32_t node_id, char* names, const uint32_t name_stride, const char* prefix, const char*suffix);
    void print_node_flattened(Node* node_data, uint32_t node_id, char* names, const uint32_t name_stride, const char* prefix, const char*suffix);

//    Node* make_node(uint32_t, Node* left, Node* right, double dist);
#ifdef __cplusplus
}
#endif