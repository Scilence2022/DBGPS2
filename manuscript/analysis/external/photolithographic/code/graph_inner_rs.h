#ifndef ANTK_GRAPH_INNER_RS_H
#define ANTK_GRAPH_INNER_RS_H
#ifdef __cplusplus
extern "C" {
#endif
void inner_begin(unsigned n, const char *prefix);
int inner_decode(const char *input, char *output, unsigned *index, int *position);
int inner_add(const char *input, unsigned source_index);
void inner_finish(const char *reference);
#ifdef __cplusplus
}
#endif
#endif
