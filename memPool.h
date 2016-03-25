/*
 * memPool.h
 *
 */

#ifndef MEMPOOL_H_
#define MEMPOOL_H_

#include <cstdlib>
#include <cmath>

class node
{
public:
	char* buf;
	node* next;
	node* prev;
private:
};

template <class T>
class MemPool
{
public:
	MemPool(unsigned long int blocksize = 0);
	~MemPool ();
	void formatReserve(size_t nelem, T initialVal = 0);
	void rawReserve (size_t byte_size);
	void deleteReserve();
	void addRaw(size_t byte_size);
	void addFormated(size_t nelem, T initialVal);
	size_t blockSize () const;
	size_t nBlocks () const;
	node* headNode () const;
	node* lastNode () const;
	int initial () const;
	int status () const;
private:
	// private variables
	mutable int fail;
	int initialized;
	size_t _nbytes; // total number of bytes to reserve
	size_t _blocksz; // max number bytes per block
	size_t _nblock; // number of blocks
	size_t _perBlockN; // number of elements per block
	node* _head; // beginning node
	node* _back; // ending node
	//private functions
	void setBlockSize (unsigned long int sz);
};

template <class T> MemPool<T>::MemPool (size_t blocksize)
	: fail(0),
	  initialized(0),
	  _nbytes(0),
	  _blocksz(10000000),
	  _nblock(0),
	  _perBlockN(0),
	  _head(0),
	  _back(0)
{
	if (blocksize != 0)
		setBlockSize(blocksize);
}

template <class T> MemPool<T>::~MemPool ()
{
	deleteReserve();
}

template <class T> void MemPool<T>::deleteReserve ()
{
	size_t sz = sizeof(T);
	node* curr_node = _head;
	node* prev_node = 0;
	while (curr_node)
	{
		if (initialized)
		{
			T* elem = 0;
			for (size_t j = 0; j < _perBlockN; ++j)
			{
				elem = reinterpret_cast <T*> (curr_node->buf + sz * j);
				elem->~T();
			}
		}
		prev_node = curr_node;
		curr_node = curr_node->next;
		delete [] prev_node->buf;
		delete prev_node;
	}
	_head = 0;
	_back = 0;
	initialized = 0;
	_nbytes = 0;
	_nblock = 0;
	_perBlockN = 0;
}

template <class T> void MemPool<T>::setBlockSize (size_t sz)
{
	_blocksz = sz;
}

template <class T> void MemPool<T>::rawReserve (size_t byte_size)
{
	_nbytes = byte_size;
	_nblock = ceil(static_cast<float>(_nbytes)/_blocksz);
	_head = new node;
	_head->prev = 0;
	node* curr_node = _head;
	node* prev_node = 0;
	for (size_t i = 0; i < _nblock; ++i)
	{
		curr_node->buf = new char[_blocksz];
		if (i == _nblock - 1)
		{
			curr_node->next = 0;
			curr_node->prev = prev_node;
			_back = curr_node;
		}
		else
		{
			prev_node = curr_node;
			curr_node->next = new node;
			curr_node = curr_node->next;
			curr_node->prev = prev_node;
		}
	}
}


template <class T> void MemPool<T>::formatReserve(size_t nelem, T initialVal)
{
	size_t sz = sizeof(T);
	_nbytes = nelem * sz;
	while (_blocksz % sz)
		++_blocksz;
	_nblock = ceil(static_cast<float>(_nbytes)/_blocksz);
	_perBlockN = _blocksz/sz;
	T* new_elem = 0;
	size_t j = 0;
	_head = new node;
	_head->prev = 0;
	node* curr_node = _head;
	node* prev_node = 0;
	for (size_t i = 0; i < _nblock; ++i)
	{
		curr_node->buf = new char[_blocksz];
		for (j = 0; j < _perBlockN; ++j)
		{
			new_elem = new (curr_node->buf + sz * j) T; // new (_head->buf + sz * j) T
			*new_elem = initialVal;
		}
		if (i == _nblock - 1 )
		{
			curr_node->next = 0;
			curr_node->prev = prev_node;
			_back = curr_node;
		}
		else
		{
			prev_node = curr_node;
			curr_node->next = new node;
			curr_node = curr_node->next;
			curr_node->prev = prev_node;
		}
	}
	initialized = 1;
}

template <class T> void MemPool<T>::addRaw(size_t byte_size)
{
	size_t newblock = ceil(byte_size/_blocksz);
	_nblock += newblock;
	node* curr_node = 0;
	node* prev_node = _back;
	prev_node->next = curr_node;

	for (size_t i = 0; i < newblock; ++i)
	{
		curr_node->buf = new char[_blocksz];
		if (i == _nblock - 1)
		{
			curr_node->next = 0;
			curr_node->prev = prev_node;
			_back = curr_node;
		}
		else
		{
			prev_node = curr_node;
			curr_node->next = new node;
			curr_node = curr_node->next;
			curr_node->prev = prev_node;
		}
	}

}

template <class T> void MemPool<T>::addFormated(size_t nelem, T initialVal)
{
	size_t sz = sizeof(T);
	size_t newbytes = nelem * sz;
	_nbytes += newbytes;
	size_t newblock = ceil(newbytes/_blocksz);
	_nblock += newblock;
	node* curr_node = 0;
	node* prev_node = _back;
	prev_node->next = curr_node;
	T* new_elem = 0;
	size_t j = 0;
	for (size_t i = 0; i < newblock; ++i)
	{
		curr_node->buf = new char[_blocksz];
		for (j = 0; j < _perBlockN; ++j)
		{
			new_elem = new (curr_node->buf + sz * j) T;
			*new_elem = initialVal;
		}
		if (i == _nblock - 1 )
		{
			curr_node->next = 0;
			curr_node->prev = prev_node;
			_back = curr_node;
		}
		else
		{
			prev_node = curr_node;
			curr_node->next = new node;
			curr_node = curr_node->next;
			curr_node->prev = prev_node;
		}
	}
}

template <class T> size_t MemPool<T>::blockSize () const
{
	return _perBlockN;
}

template <class T> size_t MemPool<T>::nBlocks () const
{
	return _nblock;
}

template <class T> node* MemPool<T>::headNode () const
{
	return _head;
}

template <class T> node* MemPool<T>::lastNode () const
{
	return _back;
}
template <class T> int MemPool<T>::initial () const
{
	return initialized;
}
template <class T> int MemPool<T>::status () const
{
	return fail;
}

#endif /* MEMPOOL_H_ */
