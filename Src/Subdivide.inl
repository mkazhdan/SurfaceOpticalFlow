/*
Copyright (c) 2016, Michael Kazhdan and Fabian Prada
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef SUBDIVIDE
#define SUBDIVIDE

template< class Real , class Vertex > int Subdivide( std::vector< TriangleIndex >& triangles , std::vector< Vertex >& vertices , Real edgeLength );
template< class Real , class Vertex , class Data > int Subdivide( std::vector< TriangleIndexWithData< Data > >& triangles , std::vector< Vertex >& vertices , Real edgeLength );

/////////////////////////////////////////////////////



template< class Real , class Vertex >
int _Subdivide( std::vector< TriangleIndex >& triangles , std::vector< Vertex >& vertices , Real edgeLength )
{
	int subdivideCount = 0;
	std::unordered_map< long long , int > edgeMap;
	std::vector< TriangleIndex > _triangles;
	std::vector< Vertex > _vertices = vertices;

	for( int i=0 ; i<triangles.size() ; i++ )
	{
		int eCount = 0;
		int e[] = { -1 , -1 , -1 };
		for( int j=0 ; j<3 ; j++ ) 
		{
			int i1 = triangles[i][j] , i2 = triangles[i][(j+1)%3];
			Real l2 = Point3D< Real >::SquareNorm( Point3D< Real >( vertices[i1] ) - Point3D< Real >( vertices[i2] ) );
			if( l2>edgeLength*edgeLength )
			{
				long long key = EdgeKey( i1 , i2 );
				int idx;
				std::unordered_map< long long , int >::iterator iter = edgeMap.find( key );
				if( iter==edgeMap.end() )
				{
					idx = (int)_vertices.size();
					edgeMap[key] = idx;
					_vertices.push_back( ( vertices[i1] + vertices[i2] ) / 2 );
					subdivideCount++;
				}
				else idx = iter->second;
				e[j] = idx;
				eCount++;
			}
		}
		if( eCount==0 ) _triangles.push_back( triangles[i] );
		else if( eCount==1 )
		{
			for( int j=0 ; j<3 ; j++ ) if( e[j]!=-1 )
			{
				_triangles.push_back( TriangleIndex( triangles[i][ j     ] , e[j] , triangles[i][(j+2)%3] ) );
				_triangles.push_back( TriangleIndex( triangles[i][(j+1)%3] , triangles[i][(j+2)%3] , e[j] ) );
			}
		}
		else if( eCount==2 )
		{
			for( int j=0 ; j<3 ; j++ ) if( e[j]==-1 )
			{
				_triangles.push_back( TriangleIndex( e[(j+1)%3] , triangles[i][(j+2)%3] , e[(j+2)%3] ) );
				_triangles.push_back( TriangleIndex( triangles[i][j] , triangles[i][(j+1)%3] , e[(j+2)%3] ) );
				_triangles.push_back( TriangleIndex( triangles[i][(j+1)%3] , e[(j+1)%3] , e[(j+2)%3] ) );
			}
		}
		else if( eCount==3 )
		{
			for( int j=0 ; j<3 ; j++ ) _triangles.push_back( TriangleIndex( triangles[i][j] , e[j] , e[(j+2)%3] ) );
			_triangles.push_back( TriangleIndex( e[0] , e[1] , e[2] ) );
		}
	}
	triangles = _triangles;
	vertices = _vertices;
	return subdivideCount;
}
template< class Real , class Vertex >
int Subdivide( std::vector< TriangleIndex >& triangles , std::vector< Vertex >& vertices , Real edgeLength )
{
	int count = 0;
	while( true )
	{
		int _count = _Subdivide( triangles , vertices , edgeLength );
		if( _count ) count += _count;
		else return count;
	}
}
template< class Real , class Vertex , class Data >
int _Subdivide( std::vector< TriangleIndexWithData< Data > >& triangles , std::vector< Vertex >& vertices , Real edgeLength )
{
	int subdivideCount = 0;
	std::unordered_map< long long , int > edgeMap;
	std::vector< TriangleIndexWithData< Data > > _triangles;
	std::vector< Vertex > _vertices = vertices;

	for( int i=0 ; i<triangles.size() ; i++ )
	{
		int eCount = 0;
		int e[] = { -1 , -1 , -1 };
		Data data[3];
		for( int j=0 ; j<3 ; j++ ) 
		{
			int i1 = triangles[i][j] , i2 = triangles[i][(j+1)%3];
			Real l2 = Point3D< Real >::SquareNorm( Point3D< Real >( vertices[i1] ) - Point3D< Real >( vertices[i2] ) );
			if( l2>edgeLength*edgeLength )
			{
				long long key = EdgeKey( i1 , i2 );
				int idx;
				std::unordered_map< long long , int >::iterator iter = edgeMap.find( key );
				if( iter==edgeMap.end() )
				{
					idx = (int)_vertices.size();
					edgeMap[key] = idx;
					_vertices.push_back( ( vertices[i1] + vertices[i2] ) / 2 );
					subdivideCount++;
				}
				else idx = iter->second;
				e[j] = idx;
				data[j] = ( triangles[i].data[j] + triangles[i].data[(j+1)%3] ) / 2;
				eCount++;
			}
		}
		if( eCount==0 ) _triangles.push_back( triangles[i] );
		else if( eCount==1 )
		{
			for( int j=0 ; j<3 ; j++ ) if( e[j]!=-1 )
			{
				_triangles.push_back( TriangleIndexWithData< Data >( triangles[i][ j     ] , e[j] , triangles[i][(j+2)%3] , triangles[i].data[ j     ] , data[j] , triangles[i].data[(j+2)%3] ) );
				_triangles.push_back( TriangleIndexWithData< Data >( triangles[i][(j+1)%3] , triangles[i][(j+2)%3] , e[j] , triangles[i].data[(j+1)%3] , triangles[i].data[(j+2)%3] , data[j] ) );
			}
		}
		else if( eCount==2 )
		{
			for( int j=0 ; j<3 ; j++ ) if( e[j]==-1 )
			{
				_triangles.push_back( TriangleIndexWithData< Data >( e[(j+1)%3] , triangles[i][(j+2)%3] , e[(j+2)%3] , data[(j+1)%3] , triangles[i].data[(j+2)%3] , data[(j+2)%3] ) );
				_triangles.push_back( TriangleIndexWithData< Data >( triangles[i][j] , triangles[i][(j+1)%3] , e[(j+2)%3] , triangles[i].data[j] , triangles[i].data[(j+1)%3] , data[(j+2)%3] ) );
				_triangles.push_back( TriangleIndexWithData< Data >( triangles[i][(j+1)%3] , e[(j+1)%3] , e[(j+2)%3] , triangles[i].data[(j+1)%3] , data[(j+1)%3] , data[(j+2)%3] ) );
			}
		}
		else if( eCount==3 )
		{
			for( int j=0 ; j<3 ; j++ ) _triangles.push_back( TriangleIndexWithData< Data >( triangles[i][j] , e[j] , e[(j+2)%3] , triangles[i].data[j] , data[j] , data[(j+2)%3] ) );
			_triangles.push_back( TriangleIndexWithData< Data >( e[0] , e[1] , e[2] , data[0] , data[1] , data[2] ) );
		}
	}
	triangles = _triangles;
	vertices = _vertices;
	return subdivideCount;
}

template< class Real , class Vertex , class Data >
int Subdivide( std::vector< TriangleIndexWithData< Data > >& triangles , std::vector< Vertex >& vertices , Real edgeLength )
{
	int count = 0;
	while( true )
	{
		int _count = _Subdivide( triangles , vertices , edgeLength );
		if( _count ) count += _count;
		else return count;
	}
}
#endif // SUBDIVIDE
