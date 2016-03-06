#ifndef MULTIVECTOR_H
#define MULTIVECTOR_H

#include "KBlade.h"

typedef struct Multivector
{
	KBlade* terms;
	unsigned int dimension;
} Multivector;

///
//Allocates a multivector
//
//Returns:
//	A pointer to a newly allocated uninitialized multivector
Multivector* Multivector_Allocate(void);

///
//Initializes a multivector
//
//Parameters:
//	mVector: A pointer to the multivector to initialize
//	dimension: the dimension of the space in which this multivector exists
void Multivector_Initialize(Multivector* mVector, int dimension);

///
//Frees memory allocated for a multivector
//
//Parameters:
//	mVector: A pointer to the multivector to free
void Multivector_Free(Multivector* mVector);

///
//Computes the sum of two multivectors
//
//Parameters:
//	dest: a pointer to the multivector to store the sum
//	a: a pointer to the LHS operand
//	b: a pointer to the RHS operand
void Multivector_Add(Multivector* dest, Multivector* a, Multivector* b);

///
//Computes the difference of two multivectors
//That is, a - b.
//
//Parameters:
//	dest: a pointer to the multivector to store the difference
//	a: a pointer to the LHS operand
//	b: a pointer to the RHS operand
void Multivector_Subtract(Multivector* dest, Multivector* a, Multivector* b);

///
//Computes the product of two multivectors,
//that is, a * b
//
//Parameters:
//	dest: A pointer to the multivector to store the product
//	a: A pointer to the LHS operand
//	b: A pointer to the RHS operand
void Multivector_GetProduct(Multivector* dest, Multivector* a, Multivector* b);

///
//Prints the contents of a multivector
//
//Parameters:
//	mVector: A pointer to the multivector to print
void Multivector_Print(Multivector* mVector);

#endif
