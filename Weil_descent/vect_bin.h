#ifndef _VECT_BIN_H_
#define _VECT_BIN_H_
#include <inttypes.h>
#include "constants.h"

/// type of boolean
typedef char _bool_t;
#define _true 1
#define _false 0


/// type of one cell of the binary vector
typedef uint64_t _vect_bin_t;

/// Size - in bits - of the binary vector
extern unsigned int _vect_bin_size;

/// Size - computed - of needed _vect_bin_t cells to store the binary
/// vector and a pointer to the next cell.
//_vect_bin_t v[_vect_bin_array_size];*/
extern unsigned int _vect_bin_array_size;

#define init(v) \
    _vect_bin_t v[__ARRAY_SIZE__]; \
    vect_bin_t_reset(v);



/// ----------------------------------- prototypes

void print_vect(_vect_bin_t *);
void print_vect_bin(_vect_bin_t *v, int pack_size);
_bool_t vect_bin_t_is_zero(_vect_bin_t *_v);
_vect_bin_t *vect_bin_t_reset(_vect_bin_t *);
_bool_t vect_bin_get_bit(_vect_bin_t *, int);
_bool_t vect_bin_flip_bit(_vect_bin_t *t, int rank);
void vect_bin_set_1(_vect_bin_t *, int);
void vect_bin_set_0(_vect_bin_t *, int);
void vect_bin_get_vect_bin_from(_vect_bin_t *, int, int, _vect_bin_t *);
_vect_bin_t *vect_bin_or(_vect_bin_t *_v1, _vect_bin_t *_v2);
_vect_bin_t *vect_bin_and(_vect_bin_t *_v1, _vect_bin_t *_v2);
_vect_bin_t *vect_bin_xor(_vect_bin_t *_v1, _vect_bin_t *_v2);
void set_vect_size(unsigned int new_size);
_vect_bin_t *vect_bin_left_shift(_vect_bin_t *_v, int bits);
_vect_bin_t *vect_bin_right_shift(_vect_bin_t *_v, int bits);
void vect_memcpy(_vect_bin_t *dest, _vect_bin_t *src, int pack_size, int offset);

#endif
