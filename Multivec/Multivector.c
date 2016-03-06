#include "Multivector.h"

#include <stdlib.h>
#include <stdio.h>


///
//Allocates a multivector
//
//Returns:
//	A pointer to a newly allocated uninitialized multivector
Multivector* Multivector_Allocate(void)
{
	return malloc(sizeof(Multivector));
}

///
//Initializes a multivector
//
//Parameters:
//	mVector: A pointer to the multivector to initialize
//	dimension: the dimension of the space in which this multivector exists
void Multivector_Initialize(Multivector* mVector, int dimension)
{
	mVector->dimension = dimension;
	mVector->terms = KBlade_AllocateArray(mVector->dimension + 1);

	for(unsigned int i = 0; i <= mVector->dimension; i++)
	{
		KBlade_Initialize(mVector->terms + i, i, mVector->dimension);
	}
}

///
//Frees memory allocated for a multivector
//
//Parameters:
//	mVector: A pointer to the multivector to free
void Multivector_Free(Multivector* mVector)
{
	KBlade_FreeArray(mVector->terms, mVector->dimension);
	free(mVector);
}

///
//Computes the sum of two multivectors
//
//Parameters:
//	dest: a pointer to the multivector to store the sum
//	a: a pointer to the LHS operand
//	b: a pointer to the RHS operand
void Multivector_Add(Multivector* dest, Multivector* a, Multivector* b)
{
	if(a->dimension == b->dimension)
	{
		if(a->dimension == dest->dimension)
		{
			for(unsigned int i = 0; i <= dest->dimension; i++)
			{
				KBlade_Add(dest->terms + i, a->terms + i, b->terms + i);
			}
		}
		else
		{
			printf("Multivector_Add failed! Operands are of different dimension. Sum not computed.\n");
		}
	}
	else
	{
		printf("Multivector_Add failed! Operands are of different dimension. Sum not computed.\n");
	}
}

///
//Computes the difference of two multivectors
//That is, a - b.
//
//Parameters:
//	dest: a pointer to the multivector to store the difference
//	a: a pointer to the LHS operand
//	b: a pointer to the RHS operand
void Multivector_Subtract(Multivector* dest, Multivector* a, Multivector* b)
{
	if(a->dimension == b->dimension)
	{
		if(a->dimension == dest->dimension)
		{
			for(unsigned int i = 0; i <= dest->dimension; i++)
			{
				KBlade_Subtract(dest->terms + i, a->terms + i, b->terms + i);
			}
		}
		else
		{
			printf("Multivector_Subtract failed! Operands are of different dimension. Difference not computed.\n");
		}
	}
	else
	{
		printf("Multivector_Subtract failed! Operands are of different dimension. difference not computed.\n");
	}
}

///
//Computes the product of two multivectors,
//that is, a * b
//
//Parameters:
//	dest: A pointer to the multivector to store the product
//	a: A pointer to the LHS operand
//	b: A pointer to the RHS operand
void Multivector_GetProduct(Multivector* dest, Multivector* a, Multivector* b)
{
	if(a->dimension == b->dimension)
	{
		if(dest->dimension == a->dimension)
		{
			for(unsigned int i = 0; i <= a->dimension; i++)
			{
				for(unsigned int j = 0; j <= b->dimension; i++)
				{
					KBlade_GetProduct(dest, a->terms + i, b->terms + j);
				}
			}
		}
		else
		{
			printf("Multivector_GetProduct failed! Destination not of proper dimension. Product not computed.");
		}
	}
	else
	{
		printf("Multivector_GetProduct failed! Operands are not of equal dimension. Product Not Computed.");
	}
}



///
//Prints the contents of a multivector
//
//Parameters:
//	mVector: A pointer to the multivector to print
void Multivector_Print(Multivector* mVector)
{
	for(unsigned int i = 0; i <= mVector->dimension; i++)
	{
		KBlade_Print(mVector->terms + i);
	}
}
