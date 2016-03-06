#include "KBlade.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "Multivector.h"
#include "Compute.h"

///
//Static Declarations
///

///
//Sets the values of the basis KBlades making up a linear combination of the KBlade with the components
//
//Parameters:
//	blade: The blade to generate it's basisKBlades
static void KBlade_InitializeBasisKBlades(KBlade* blade);

///
//Recursive helper function for generating basis KBlades
//Determines the bit offsets for the next KBlade to be generated
//
//Parameters:
//	e: An array of the bit offsets in the basis K Blade being generated
//	grade: The grade of the KBlade being initialized
//	dim: The dimension of the KBlade being initialized
//
//Returns:
//	A character indicating whether or not any more KBlades can be generated
static unsigned char KBlade_InitializeBasisKBladeCarry(int* e, unsigned int grade, unsigned int dim);

///
//Determines if the product of two basis K blades will negate the sign
//of the coefficient of the product
//
//Parameters:
//	basisA: A bitmask indicating the basis K blade for the LHS operand
//	gradeA: The grade of the basis K blade of the LHS operand
//	basisB: A bitmask indicating the basis K blade for the RHS operand
//	gradeB: The grade of the basisKBlade of the RHS operand
//	dimension: The dimension of the two basis KBlades
static float KBlade_DetermineProductTermSign(char basisA, unsigned int gradeA, char basisB, unsigned int gradeB, unsigned int dimension);

///
//Determines the index of a basis K-Blade term in contiguous memory
//
//Parameters:
//	blade: bit encoding of a basis KBlade
//	grade: the grade of the blade
//	dim: The dimension of the blade
static unsigned int KBlade_DetermineBasisKBladeIndex(char blade, unsigned int grade, unsigned int dim);

///
//Determines the grade of a Basis KBlade given it's bitmask
//
//Parameters:
//	v: The bitmask of the Basis KBlade
static unsigned char KBlade_GetBasisKBladeGrade(unsigned char v);

///
//Allocates a KBlade
//
//returns:
//	A pointer to a newly allocated uninitialized KBlade
KBlade* KBlade_Allocate(void)
{
	return malloc(sizeof(KBlade));
}

///
//Allocates an array of KBlades
//
//Parameters:
//	size: The size of the array to allocate
//
//Returns:
//	A pointer to an array of KBlades of length size
KBlade* KBlade_AllocateArray(int size)
{
	return malloc(sizeof(KBlade) * size);
}

///
//Initializes a KBlade
//
//Parameters:
//	blade: A pointer to the blade to initialize
//	grade: The grade of the blade
//	dimension: The dimension in which this blade lives
void KBlade_Initialize(KBlade* blade, int grade, int dimension)
{
	blade->dimension = dimension;
	blade->grade = grade;

	blade->numComponents = Compute_BinomialCoefficient(dimension, grade);
	if(blade->numComponents == 0)
	{
		blade->numComponents = 1;
	}
	blade->components = calloc(sizeof(float), blade->numComponents);

	blade->basisKBlades = calloc(sizeof(char), blade->numComponents);
	KBlade_InitializeBasisKBlades(blade);
}

///
//Free memory allocated for a KBlade
//
//Parameters:
//	blade: A pointer to the blade being freed
void KBlade_Free(KBlade* blade)
{
	free(blade->components);
	free(blade->basisKBlades);
	free(blade);
}

///
//Free memory allocated for an array of KBlades
//
//Parameters:
//	blades: A pointer to the array of blades to free
//	size: The number of blades in the array
void KBlade_FreeArray(KBlade* blades, int size)
{
	for(int i = 0; i < size; i++)
	{
		free(blades[i].components);
		free(blades[i].basisKBlades);
	}
	free(blades);
}

///
//Recursive helper function for generating basis KBlades
//Determines the bit offsets for the next KBlade to be generated
//
//Parameters:
//	e: An array of the bit offsets in the basis K Blade being generated
//	grade: The grade of the KBlade being initialized
//	dim: The dimension of the KBlade being initialized
//
//Returns:
//	A character indicating whether or not any more KBlades can be generated
static unsigned char KBlade_InitializeBasisKBladeCarry(int* e, unsigned int grade, unsigned int dim)
{
	//Confirm carry is needed
	int sGrade = grade;
	int sDim = dim;
	if(sGrade - 2 < 0) return 0;
	if(++e[sGrade - 2] >= sDim-1) 
	{
		if(grade - 1 == 1) return 0;
		return KBlade_InitializeBasisKBladeCarry(e, grade - 1, dim-1);
	}
	e[grade - 1] = e[grade - 2] + 1;
	return 1;
}

///
//Sets the values of the basis KBlades making up a linear combination of the KBlade with the components
//
//Parameters:
//	blade: The blade to generate it's basisKBlades
static void KBlade_InitializeBasisKBlades(KBlade* blade)
{
	if(blade->grade == 0)
		return;
	int* e = calloc(sizeof(int), blade->grade);
	for(unsigned int i = 0; i < blade->grade; i++)
	{
		e[i] = i;
	}

	//e = [0, 1, 2, .., Grade - 1]
	unsigned char keepGoing = 1;
	int kBladesGenerated = 0;
	while(keepGoing)
	{
		for(unsigned char i = 0; i < blade->grade; i++)
		{
			blade->basisKBlades[kBladesGenerated] |= (1 << e[i]);
		}

		kBladesGenerated++;
		
		if(++e[blade->grade - 1] >= (int)blade->dimension)
		{
			keepGoing = KBlade_InitializeBasisKBladeCarry(e, blade->grade, blade->dimension);
		}
	}

	free(e);
}


///
//Alters the grade of a KBlade and resizes the blade accordingly
//The blade will be left as the zero KBlade with that grade
//
//Parameters:
//	blade: A pointer to the KBlade to regrade
//	grade: The new grade of the KBlade
void KBlade_ReGrade(KBlade* blade, int grade)
{
	//Determine the new number of components
	int numComps = Compute_BinomialCoefficient(blade->dimension, grade);
	blade->grade = grade;

	if(numComps != blade->numComponents)
	{
		free(blade->components);
		free(blade->basisKBlades);

		blade->numComponents = numComps;
		blade->components = calloc(sizeof(float), blade->numComponents);
		blade->basisKBlades = calloc(sizeof(char), blade->numComponents);
	}
	else
	{
		memset(blade->components, 0, blade->numComponents * sizeof(float));
		memset(blade->basisKBlades, 0, blade->numComponents * sizeof(char));
	}

	KBlade_InitializeBasisKBlades(blade);

}

///
//Increments one KBlade by another
//
//Parameters:
//	dest: A pointer to the KBlade to increment
//	inc: A pointer to the KBlade to increment it by
void KBlade_Increment(KBlade* dest, KBlade* inc)
{
	if(dest->grade == inc->grade && dest->dimension == inc->dimension)
	{
		for(int i = 0; i < dest->numComponents; i++)
		{
			dest->components[i] += inc->components[i];
		}
	}
}

///
//Computes the sum of two KBlades
//
//Parameters:
//	dest: destination pointer to store sum of a and b
//	a: Pointer to LHS operand
//	b: pointer to RHS operand
void KBlade_Add(KBlade* dest, KBlade* a, KBlade* b)
{
	if(a->grade == b->grade && a->dimension == b->dimension)
	{
		if(dest->grade == a->grade && dest->dimension == a->dimension)
		{
			for(int i = 0; i < dest->numComponents; i++)
			{
				dest->components[i] = a->components[i] + b->components[i];
			}
		}
		else
		{
			printf("KBlade_Add Failed! Destination not of proper grade or dimension. Destination was not populated.");
		}
	}
	else
	{
		printf("KBlade_Add Failed! KBlades not of equal grade or dimension. Destination was not populated.\n");
	}
}

///
//Deccrements one KBlade by another
//
//Parameters:
//	dest: A pointer to the KBlade to decrement
//	dec: A pointer to the KBlade to decrement it by
void KBlade_Decrement(KBlade* dest, KBlade* dec)
{
	if(dest->grade == dec->grade && dest->dimension == dec->dimension)
	{
		for(int i = 0; i < dest->numComponents; i++)
		{
			dest->components[i] -= dec->components[i];
		}
	}

}

///
//Computes the difference of a and b
//that is, a - b.
//
//Parameters:
//	dest: destination pointer to store the difference of a and b
//	a: Pointer to the LHS operand
//	b: Pointer to the RHS operand
void KBlade_Subtract(KBlade* dest, KBlade* a, KBlade* b)
{
	if(a->grade == b->grade && a->dimension == b->dimension)
	{
		if(dest->grade == a->grade && dest->dimension == a->dimension)
		{
			for(int i = 0; i < dest->numComponents; i++)
			{
				dest->components[i] = a->components[i] - b->components[i];
			}

		}
		else
		{
			printf("KBlade_Subtract Failed! Destination not of proper grade or dimension. Destination was not populated.");
		}
	}
	else
	{
		printf("KBlade_Subtract Failed! KBlades not of equal grade or dimension. Destination was not populated.\n");
	}
}

///
//Determines if the product of two basis K blades will negate the sign
//of the coefficient of the product
//
//Parameters:
//	basisA: A bitmask indicating the basis K blade for the LHS operand
//	gradeA: The grade of the basis K blade of the LHS operand
//	basisB: A bitmask indicating the basis K blade for the RHS operand
//	gradeB: The grade of the basisKBlade of the RHS operand
//	dimension: The dimension of the two basis KBlades
static float KBlade_DetermineProductTermSign(char basisA, unsigned int gradeA, char basisB, unsigned int gradeB, unsigned int dimension)
{
	int totalBlades = gradeA + gradeB;
	char* orderedDimensions = malloc(sizeof(char) * totalBlades);
	int index = 0;
	for(unsigned int i = 0; i < dimension; i++)
	{
		if(basisA & (1 << i)) orderedDimensions[index++] = i;
	}
	for(unsigned int i = 0; i < dimension; i++)
	{
		if(basisB & (1 << i)) orderedDimensions[index++] = i;
	}


	int power = 0;	

	unsigned int i = 0;
	unsigned int j = gradeA;
	unsigned int offset = 0;
	while((i < gradeA + offset) && (j < gradeB + gradeA))
	{
		if(orderedDimensions[i] == orderedDimensions[j])
		{
			power += (j - i) - 1;

			char save = orderedDimensions[i + 1];
			memmove(orderedDimensions + i + 1, orderedDimensions + i + 2,  (j - (i + 2)) + 1);
			orderedDimensions[j] = save;

			i+=2;
			++j;
			++offset;
		}
		else if(orderedDimensions[j] < orderedDimensions[i])
		{

			power += (j - i);

			char save = orderedDimensions[i];
			memmove(orderedDimensions + i, orderedDimensions + i + 1, j - i);
			orderedDimensions[j] = save;


			++i;
			++j;
			++offset;

		}
		else
		{
			++i;
		}
	}

	if(power > 0)
	{
		return powf(-1.0f, (float)power);
	}
	return 1.0f;
}

///
//Determines the index of a basis K-Blade term in contiguous memory
//
//Parameters:
//	blade: bit encoding of a basis KBlade
//	grade: the grade of the blade
//	dim: The dimension of the blade
static unsigned int KBlade_DetermineBasisKBladeIndex(char blade, unsigned int grade, unsigned int dim)
{
	unsigned int total = 0;
	unsigned int numSet = 0;
	for(unsigned int i = 0; i < dim; i++)
	{
		if(numSet == grade)
			break;

		if(!(blade & (1 << i)))
		{
			if(grade >= i)
				++total;
			else
				total += Compute_BinomialCoefficient(dim - i, grade - i);
		}
		else
		{
			numSet++;
		}
	}
	return total;
}

///
//Determines the grade of a Basis KBlade given it's bitmask
//
//Parameters:
//	v: The bitmask of the Basis KBlade
static unsigned char KBlade_GetBasisKBladeGrade(unsigned char v)
{
       static const unsigned char lookupTable[] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};
       return lookupTable[v & 0x0F] + lookupTable[v >> 4];
}

///
//Gets the dual of a blades
// 
//Parameters:
//	dest: A pointer to a KBlade to store the dual
//	blade: A pointer to a KBlade to take a dual of
void KBlade_GetDual(KBlade* dest, KBlade* blade)
{
	if(dest->dimension == blade->dimension)
	{
		//Calculate and check grade
		unsigned int grade = blade->dimension - blade->grade;
		if(grade == dest->grade)
		{
			//Generate the psuedo scalar for the blades dimension
			unsigned char basisPsuedoScalar = 0;
			for(unsigned int i = 0; i < blade->dimension; i++)
			{
				basisPsuedoScalar |= (1 << i);
			}

			//Multiply each term
			for(unsigned int i = 0; i < blade->dimension; i++)
			{
				float sign = KBlade_DetermineProductTermSign(blade->basisKBlades[i], blade->grade, basisPsuedoScalar, blade->dimension, blade->dimension);
				unsigned char basisKBlade = blade->basisKBlades[i] ^ basisPsuedoScalar;
				unsigned int index = KBlade_DetermineBasisKBladeIndex(basisKBlade, grade, dest->dimension);
				dest->components[index] += sign * blade->components[i];

			}

		}
		else
		{
			printf("Error: KBlade_GetDual failed! Destination is not of proper grade! Dual not computed.\n");
		}
	

	}
	else
	{
		printf("Error: KBlade_GetDual Failed! Destination is not of proper dimension. Dual not computed.\n");
	}
}

///
//Gets the geometric product of two KBlades
//
//Parameters:
//	dest: A pointer to the destination of the multivector to store the product
//	a: A pointer to the LHS operand
//	b: A pointer to the RHS operand
void KBlade_GetProduct(Multivector* dest, KBlade* a, KBlade* b)
{
	if(a->dimension == b->dimension)
	{
		if(dest->dimension == a->dimension)
		{
			for(int i = 0; i < a->numComponents; i++)
			{
				for(int j = 0; j < b->numComponents; j++)
				{
					float sign = KBlade_DetermineProductTermSign(a->basisKBlades[i], a->grade, b->basisKBlades[j], b->grade, dest->dimension);
					float scalar = a->components[i] * b->components[j];
					char basisKBlade = a->basisKBlades[i] ^ b->basisKBlades[j];
					unsigned int grade = KBlade_GetBasisKBladeGrade(basisKBlade);
					unsigned int index = KBlade_DetermineBasisKBladeIndex(basisKBlade, grade, dest->dimension);
					KBlade* dBlade = dest->terms + grade;
					dBlade->components[index] += sign * scalar;

				}
			}
		}
		else
		{
			printf("Error: KBlade_GetProduct failed! Destination is not of proper dimension! Product not computed!\n");
		}
	}
	else
	{
		printf("Error: KBlade_GetProduct failed! Operands are not of equal dimension! Product not computed!\n");
	}
}

///
//Prints the contents of a K-Blade to standard output
//
//Parameters:
//	blade: A pointer to the blade to print
void KBlade_Print(KBlade* blade)
{
	for(int i = 0; i < blade->numComponents; i++)
	{
		printf("%f e", blade->components[i]);
		for(unsigned int j = 0; j < blade->dimension; j++)
		{
			if(blade->basisKBlades[i] & (1 << j))
			{
				printf("%d", j + 1);
			}
		}

		printf(", ");
	}
	printf("\n");
}

/*
///
//Performs the inner product on two KBlades
//
//Parameters:
//	dest: a pointer to a KBlade to store the product in
//	a: a pointer to the LHS KBlade
//	b: a pointer to the RHS KBlade
void KBlade_InnerProduct(KBlade* dest, KBlade* a, KBlade* b)
{

	int maxGrade = a->grade > b->grade ? a->grade : b->grade;
	int minGrade = maxGrade == a->grade ? b->grade : a->grade;
	KBlade_ReGrade(dest, maxGrade - minGrade);

	for(int i = 0; i <

}
*/


