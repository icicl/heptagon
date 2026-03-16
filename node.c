#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef struct Node Node;
struct Node {
    uint32_t id;
    uint32_t left;
    uint32_t right;
    uint32_t parent;
    double dist;
};

Node* node_data;

void print_node(Node* node_data, uint32_t node_id, char* names, const uint32_t name_stride, const char* prefix, const char* suffix) {
    printf("%s", prefix);
    Node node = node_data[node_id];
    if (node.left == UINT32_MAX) {
        printf("%s", names + node.id*name_stride);
    } else {
        print_node(node_data, node.left, names, name_stride, "(\n", ",\n");
        print_node(node_data, node.right, names, name_stride, "", ")\n");
    }
    printf(":%.5lf", (node.dist > 0.0) ? node.dist : 0.0);
    printf("%s", suffix);
}

void print_node_flattened(Node* node_data, uint32_t node_id, char* names, const uint32_t name_stride, const char* prefix, const char* suffix) {
    printf("%s", prefix);
    Node node = node_data[node_id];
    if (node.left == UINT32_MAX) {
        printf("%s", names + node.id*name_stride);
    } else {
        print_node_flattened(node_data, node.left, names, name_stride, "(", ",");
        print_node_flattened(node_data, node.right, names, name_stride, "", ")");
    }
    printf(":%.5lf", (node.dist > 0.0) ? node.dist : 0.0);
    printf("%s", suffix);
}

/*Node* make_node(uint32_t id, Node* left, Node* right, double dist) {
    Node* node = malloc(sizeof(Node));
    node->id = id;
    node->left = left;
    node->right = right;
    node->dist = dist;
    return node;
}*/
