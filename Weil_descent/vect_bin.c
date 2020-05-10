#include<assert.h>
#include<stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <stdint.h>
#include "vect_bin.h"

/// Adressing bytes (and bits) is like this in _vect_bin_t type
///      0         1        2         3         4
/// [N-------][--------][--------][--------][-------0]
/// [--------][--------][--------][--------][--------]

const _vect_bin_t masks[64] =
{
	1ULL,2ULL,4ULL,8ULL,16ULL,32ULL,64ULL,128ULL,256ULL,512ULL,1024ULL,2048ULL,4096ULL,8192ULL,16384ULL,32768ULL,65536ULL,131072ULL,262144ULL,524288ULL,1048576ULL,2097152ULL,4194304ULL,8388608ULL,16777216ULL,33554432ULL,67108864ULL,134217728ULL,268435456ULL,536870912ULL,1073741824ULL,2147483648ULL,4294967296ULL,8589934592ULL,17179869184ULL,34359738368ULL,68719476736ULL,137438953472ULL,274877906944ULL,549755813888ULL,1099511627776ULL,2199023255552ULL,4398046511104ULL,8796093022208ULL,17592186044416ULL,35184372088832ULL,70368744177664ULL,140737488355328ULL,281474976710656ULL,562949953421312ULL,1125899906842624ULL,2251799813685248ULL,4503599627370496ULL,9007199254740992ULL,18014398509481984ULL,36028797018963968ULL,72057594037927936ULL,144115188075855872ULL,288230376151711744ULL,576460752303423488ULL,1152921504606846976ULL,2305843009213693952ULL,4611686018427387904ULL,9223372036854775808ULL
};

unsigned int _vect_bin_size;
    
unsigned int _vect_bin_array_size;

/// get bit at rank return _true or _false respectively to 1 and 0
inline _bool_t vect_bin_get_bit(_vect_bin_t *t, int rank) {
//  printf(" Rank(%d) : t[%lu] @ Rank(%lu)\n", rank, _vect_bin_array_size - 1 - (rank / (sizeof(_vect_bin_t) << 3)), (rank % (sizeof(_vect_bin_t) << 3)));
	return((t[_vect_bin_array_size - 1 - (rank / (sizeof(_vect_bin_t) << 3))] & masks[rank % (sizeof(_vect_bin_t) << 3)]) ? 1 : 0);
}

inline _bool_t vect_bin_flip_bit(_vect_bin_t *t, int rank) {
	return((t[_vect_bin_array_size - 1 - (rank / (sizeof(_vect_bin_t) << 3))] ^= masks[rank % (sizeof(_vect_bin_t) << 3)]) ? 1 : 0);
}

/// set 'rank' bit to 1
inline void vect_bin_set_1(_vect_bin_t *t, int rank) {
	t[_vect_bin_array_size - 1 - (rank / (sizeof(_vect_bin_t) << 3))] |= masks[rank % (sizeof(_vect_bin_t) << 3)];
}

/// set 'rank' bit to 0
inline void vect_bin_set_0(_vect_bin_t *t, int rank) {
  t[_vect_bin_array_size - 1 - (rank / (sizeof(_vect_bin_t) << 3))] &= (~(masks[rank % (sizeof(_vect_bin_t) << 3)]));
}

void print_vect(_vect_bin_t *v) {
  for(int i = 0; i < _vect_bin_array_size; ++i) printf(" %llu", v[i]);
  printf("\n");
}

void print_vect_bin(_vect_bin_t *v, int pack_size) {
  for(int i = 0; i < _vect_bin_size; ++i) {
	  printf("%d ", vect_bin_get_bit(v, i) ? 1 : 0);
      if(i % pack_size == pack_size - 1)
		  //printf("  ");
		  printf("\n");

  }
  printf("\n");
	printf("\n");
}

/// check if _vect_bin_t is assigned to 0
_bool_t vect_bin_t_is_zero(_vect_bin_t *_v) {
	for(int i = 0; i < _vect_bin_array_size; ++i)
		if(_v[i] != 0) return _false;
	return _true;
}

/// _vect_bin_t assigned to 0
_vect_bin_t *vect_bin_t_reset(_vect_bin_t *_v) {
  if(_v == NULL) return(NULL);
  for(int i = 0; i < _vect_bin_array_size; ++i) _v[i] = 0;
  return _v;
}

/// return a substring from t from from_bit to to_bit. o could be NULL. It's then allocated
void vect_bin_get_vect_bin_from(_vect_bin_t *t, int from_bit, int to_bit, _vect_bin_t *o) {
  int i, j;
  for(i = from_bit, j = 0; i < to_bit; ++i, ++j)
    if(vect_bin_get_bit(t, i)) vect_bin_set_1(o, j);
    else vect_bin_set_0(o, j);
}

/// v1 <- v1 OR v2
_vect_bin_t *vect_bin_or(_vect_bin_t *_v1, _vect_bin_t *_v2) {
	if(_v1 == NULL) return(NULL);
	if(_v2 == NULL) return(NULL);
	for(int i = 0; i < _vect_bin_array_size; ++i) _v1[i] |= _v2[i];
	return _v1;
}

/// v1 <- v1 AND v2
_vect_bin_t *vect_bin_and(_vect_bin_t *_v1, _vect_bin_t *_v2) {
	if(_v1 == NULL) return(NULL);
	if(_v2 == NULL) return(NULL);
	for(int i = 0; i < _vect_bin_array_size; ++i) _v1[i] &= _v2[i];
	return _v1;
}

/// v1 <- v1 XOR v2
_vect_bin_t *vect_bin_xor(_vect_bin_t *_v1, _vect_bin_t *_v2) {
  if(_v1 == NULL) return(NULL);
  if(_v2 == NULL) return(NULL);
  for(int i = 0; i < _vect_bin_array_size; ++i) _v1[i] ^= _v2[i];
  return _v1;
}

/// left shift (att: size of vector after shift should not exceed total vector size)
_vect_bin_t *vect_bin_left_shift(_vect_bin_t *_v, int bits) {
  if(_v == NULL) return(NULL);
	if(bits == 0) return _v;
    int i;
    int div = bits / (sizeof(_vect_bin_t) << 3);
    char mod = bits % (sizeof(_vect_bin_t) << 3);
    char _mod = (sizeof(_vect_bin_t) << 3) - (bits % (sizeof(_vect_bin_t) << 3));
	if(mod == 0) //because _mod==64 and shift >= width type -> undefined behavior
	{
		for (i = div; i < _vect_bin_array_size - 1; i++)
		{
			assert(i-div >= 0);
			assert(i >= 0);
			assert(i-div < _vect_bin_array_size);
			assert(i < _vect_bin_array_size);
			_v[i - div] = _v[i];
		}
	}
	else
	{
		for (i = div; i < _vect_bin_array_size - 1; i++)
		{
			assert(i-div >= 0);
			assert(i >= 0);
			assert(i-div < _vect_bin_array_size);
			assert(i + 1 >= 0);
			assert(i+ 1 < _vect_bin_array_size);
			assert(i < _vect_bin_array_size);
			_v[i - div] = (_v[i] << mod) | (_v[i + 1] >> _mod);
		}
	}
	assert(i-div >= 0);
	assert(i >= 0);
	assert(i-div < _vect_bin_array_size);
	assert(i < _vect_bin_array_size);
	_v[i - div] = (_v[i] << mod);
	div--;
    while(div >= 0)
    {
		assert(i-div >= 0);
		assert(i-div < _vect_bin_array_size);
        _v[i - div] = 0;
        div--;
    }
  return _v;
}

_vect_bin_t *vect_bin_right_shift(_vect_bin_t *_v, int bits) {
	if(_v == NULL) return(NULL);
	if(bits == 0) return _v;
	int i;
	int div = bits / (sizeof(_vect_bin_t) << 3);
	int mod = bits % (sizeof(_vect_bin_t) << 3);
	int _mod = (sizeof(_vect_bin_t) << 3) - (bits % (sizeof(_vect_bin_t) << 3));
	if(mod == 0) //because _mod==64 and shift >= width type -> undefined behavior
	{
		for (i = _vect_bin_array_size - 1 - div; i > 0 ; i--)
		{
			assert(i+div >= 0);
			assert(i >= 0);
			assert(i+div < _vect_bin_array_size);
			assert(i < _vect_bin_array_size);
			_v[i + div] = _v[i];
		}
	}
	else
	{
		for (i = _vect_bin_array_size - 1 - div; i > 0 ; i--)
		{
			//printf("%d <- v[i]-%llu %llu >> %llu = %llu.....v[i-1]-%llu %llu << %llu = %llu...=%llu\n",i + div,i,_v[i],mod,_v[i] >> mod,i-1,_v[i-1],_mod,_v[i-1] << _mod,(_v[i] >> mod) | (_v[i - 1] << _mod));
			assert(i+div >= 0);
			assert(i >= 0);
			assert(i-1 < _vect_bin_array_size);
			assert(i -1>= 0);
			assert(i+div < _vect_bin_array_size);
			assert(i < _vect_bin_array_size);
			_v[i + div] = (_v[i] >> mod) | (_v[i - 1] << _mod);
		}
	}
	assert(i+div >= 0);
	assert(i >= 0);
	assert(i+div < _vect_bin_array_size);
	assert(i < _vect_bin_array_size);
	_v[i + div] = (_v[i] >> mod);
	div--;
	while(div >= 0)
	{
		assert(i+div >= 0);
		assert(i+div < _vect_bin_array_size);
		_v[i + div] = 0;
		div--;
	}
	return _v;
}

void set_vect_size(unsigned int new_size) {
    _vect_bin_size = new_size;
    
    _vect_bin_array_size =
    ((_vect_bin_size / (sizeof(_vect_bin_t) << 3)) +
    ((_vect_bin_size % (sizeof(_vect_bin_t) << 3)) ? 1 : 0));
	
	if(_vect_bin_array_size > __ARRAY_SIZE__)
		printf("ARRAY_SIZE should be %d\n", _vect_bin_array_size);
	assert(_vect_bin_array_size <= __ARRAY_SIZE__);
}
