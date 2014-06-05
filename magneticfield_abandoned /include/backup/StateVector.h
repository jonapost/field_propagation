/*
 * StateVector.h
 *
 *  Created on: Mar 4, 2014
 *      Author: swenzel
 *
 *  Revised to adapt Geant4
 *  change int, double to G4int, G4double
 */

#ifndef STATEVECTOR_H_
#define STATEVECTOR_H_

#include "G4Types.hh"
#include "globals.hh"

// this is a physics state vector for a Runge Kutta solver
// the first  N entries describe generalized coordinates "x"
// the second N entries describe generalized impulses "q"

// this class is not optimized for anything; it should just demonstrate the concept of an abstraction
template
<G4int N>
class
State
{
public:
	// with this field, solvers or physics can infer the dimension of the problem
	static const G4int size = N;


	G4double x[N]; // the coordinates
	G4double xp[N]; // the time derivates of x

// constructor
	State(){
		for(G4int i=0;i<N;i++)
		{
			x[i]=0.;xp[i]=0.;
		}
	}

	// copy constructor ( should actually use copy and swap idiom )
	State( State<N> const & rhs ){
		for(G4int i=0;i<N;i++)
				{
					x[i]=rhs.x[i];xp[i]=rhs.xp[i];
				}
	}

	State<N> & operator+=( State<N> const & rhs )
	{
		for(G4int i=0;i<N;i++)
		{
			this->x[i]+=rhs.x[i];
			this->xp[i]+=rhs.xp[i];
		}
		return *this;
	}

	State<N> & operator=( State<N> const & rhs )
	{
		for(G4int i=0;i<N;i++)
		{
			this->x[i]=rhs.x[i];this->xp[i]=rhs.xp[i];
		}
		return *this;
	}

// operator +
	template<G4int M>
	friend
	State<M> operator+( State<M> const & lhs, State<M> const & rhs );

	template<G4int M>
	friend
	std::ostream & operator<<( std::ostream & str, State<M> const & s );

	template<G4int M>
	friend
	State<M> operator*( G4double , State<M> const & rhs );

};

template <G4int M>
std::ostream & operator<<( std::ostream & str, State<M> const & s )
{
	for(G4int i=0;i<M;i++)
			{
				str << s.x[i] << " ";
			}
	for(G4int i=0;i<M;i++)
			{
				str << s.xp[i] << " ";
			}
	return str;
}


template<G4int N>
inline
State<N> operator+( State<N> const &lhs, State<N> const & rhs )
{
	State<N> tmp( lhs );
	tmp+=rhs;
	return tmp;
}

template<G4int N>
inline
State<N> operator*( G4double scalar, State<N> const & rhs )
{
	State<N> tmp( rhs );
	for(G4int i=0;i<N;i++)
	{
		tmp.x[i]*=scalar;
		tmp.xp[i]*=scalar;
	}
	return tmp;
}

#endif /* STATEVECTOR_H_ */
