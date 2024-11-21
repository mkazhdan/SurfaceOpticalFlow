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

#ifndef TEXTURE_INCLUDED
#define TEXTURE_INCLUDED
#include <Misha/FEM.h>

template< class Real >
typename FEM::RiemannianMesh< Real >::SamplePoint* GetTextureSource
(
	const FEM::RiemannianMesh< Real >* mesh ,
	ConstPointer( typename FEM::RiemannianMesh< Real >::CoordinateXForm ) xForms ,
	const std::vector< Point2D< Real > >& triangleTextures ,
	int tWidth , int tHeight ,
	int padRadius
);

template< class Real , class Vertex >
void SampleTextureToVertices
(
	const std::vector< TriangleIndexWithData< Point2D< Real > > >& triangles ,
	const std::vector< Vertex >& vertices , 
	const unsigned char* texture , int tWidth , int tHeight , 
	std::vector< Point3D< Real > >& colors ,
	bool antiAliased=true
);

template< class Real >
Point3D< Real > SampleTexture
(
	const unsigned char* texture ,
	int tWidth , int tHeight ,
	Point2D< Real > p ,
	bool bilinear=true
);

template< class Real >
void RasterizeTriangleInterior
(
	Point2D< Real > v0 , Point2D< Real > v1 , Point2D< Real > v2 ,
	int tIdx ,
	FEM::SamplePoint< Real >* samplePoints ,
	int width , int height
);

/////////////////////////////////////////////////////

template< class Real >
Point3D< Real > SampleTexture( const unsigned char* texture , int tWidth , int tHeight , Point2D< Real > p , bool bilinear )
{
	p[1] = 1 - p[1];
	p[0] = std::min< Real >( 1. , std::max< Real >( 0. , p[0] ) );
	p[1] = std::min< Real >( 1. , std::max< Real >( 0. , p[1] ) );
	p[0] *= tWidth-1 , p[1] *= tHeight-1;
	int x0 = (int)floor( p[0] ) , y0 = (int)floor( p[1] );
	if( bilinear )
	{
		Real dx = p[0] - x0 , dy = p[1] - y0;
		int x1 = std::min< int >( x0+1 , tWidth-1 ) , y1 = std::min< int >( y0+1 , tHeight-1 );
		return
			Point3D< Real >( (Real)( texture[3*(tWidth*y0+x0)] ) , (Real)( texture[3*(tWidth*y0+x0)+1] ) , (Real)( texture[3*(tWidth*y0+x0)+2] ) ) * (Real)( (1.-dx) * (1.-dy) ) +
			Point3D< Real >( (Real)( texture[3*(tWidth*y0+x1)] ) , (Real)( texture[3*(tWidth*y0+x1)+1] ) , (Real)( texture[3*(tWidth*y0+x1)+2] ) ) * (Real)( (   dx) * (1.-dy) ) +
			Point3D< Real >( (Real)( texture[3*(tWidth*y1+x1)] ) , (Real)( texture[3*(tWidth*y1+x1)+1] ) , (Real)( texture[3*(tWidth*y1+x1)+2] ) ) * (Real)( (   dx) * (   dy) ) +
			Point3D< Real >( (Real)( texture[3*(tWidth*y1+x0)] ) , (Real)( texture[3*(tWidth*y1+x0)+1] ) , (Real)( texture[3*(tWidth*y1+x0)+2] ) ) * (Real)( (1.-dx) * (   dy) ) ;
	}
	else return Point3D< Real >( (Real)( texture[3*(tWidth*y0+x0)] ) , (Real)( texture[3*(tWidth*y0+x0)+1] ) , (Real)( texture[3*(tWidth*y0+x0)+2] ) );
}

template< class Real , class Vertex >
void SampleTextureToVertices
(
	const std::vector< TriangleIndexWithData< Point2D< Real > > >& triangles ,
	const std::vector< Vertex >& vertices , 
	const unsigned char* texture , int tWidth , int tHeight , 
	std::vector< Point3D< Real > >& colors ,
	bool antiAliased
)
{
	colors.resize( vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) colors[i] = Point3D< Real >();
	if( antiAliased )
	{
		std::vector< Real > count( vertices.size() , 0 );
		std::vector< Real > areaRatios( triangles.size() , (Real)0. );
#pragma omp parallel for
		for( int i=0 ; i<triangles.size() ; i++ )
		{
			Point3D< Real > v[] = { ( Point3D< Real > )vertices[ triangles[i][0] ] , ( Point3D< Real > )vertices[ triangles[i][1] ] , ( Point3D< Real > )vertices[ triangles[i][2] ] };
			Point3D< Real > w[3];
			for( int j=0 ; j<3 ; j++ ) w[j] = Point3D< Real >( triangles[i].data[j][0] , triangles[i].data[j][1] , (Real)0 );
			areaRatios[i] = (Real)sqrt( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ).squareNorm() / Point3D< Real >::CrossProduct( w[1]-w[0] , w[2]-w[0] ).squareNorm() );
		}

		FEM::SamplePoint< Real >* samplePoints = new FEM::SamplePoint< Real >[ tWidth * tHeight ];
		for( int i=0 ; i<tWidth*tHeight ; i++ ) samplePoints[i].tIdx = -1;
		for( int i=0 ; i<triangles.size() ; i++ )
		{
			Point2D< Real > d[] = { triangles[i].data[0] , triangles[i].data[1] , triangles[i].data[2] };
			for( int j=0 ; j<3 ; j++ ) d[j][1] = (Real)1.-d[j][1];
			RasterizeTriangleInterior( d[0] , d[1] , d[2] , i , samplePoints , tWidth , tHeight );
		}
		for( int i=0 ; i<tWidth*tHeight ; i++ ) if( samplePoints[i].tIdx!=-1 )
		{
			int j;
			Point3D< Real > b( (Real)1.-samplePoints[i].p[0]-samplePoints[i].p[1] , samplePoints[i].p[0] , samplePoints[i].p[1] );
			if     ( b[0]>=b[1] && b[0]>=b[2] ) j = 0;
			else if( b[1]>=b[0] && b[1]>=b[2] ) j = 1;
			else                                j = 2;
			colors[ triangles[ samplePoints[i].tIdx ][j] ] += Point3D< Real >( (Real)texture[3*i+0] , (Real)texture[3*i+1] , (Real)texture[3*i+2] ) * areaRatios[ samplePoints[i].tIdx ];
			count [ triangles[ samplePoints[i].tIdx ][j] ] += areaRatios[ samplePoints[i].tIdx ];
		}
		delete[] samplePoints;

		// If the vertex is surrounded by triangles whose interiors are not rasterized, default to sampling.
		{
			std::vector< Point3D< Real > > _colors( vertices.size() , Point3D< Real >() );
			std::vector< int > _count( vertices.size() , 0 );
			for( int i=0 ; i<triangles.size() ; i++ ) for( int j=0 ; j<3 ; j++ ) _colors[ triangles[i][j] ] += SampleTexture( texture , tWidth , tHeight , triangles[i].data[j] , true ) , _count[ triangles[i][j] ]++;
			for( int i=0 ; i<vertices.size() ; i++ )
				if( count[i] ) colors[i] /= (Real)count[i];
				else colors[i] = _colors[i] / (Real)_count[i];
		}
	}
	else
	{
		std::vector< int > count( vertices.size() , 0 );
		for( int i=0 ; i<triangles.size() ; i++ ) for( int j=0 ; j<3 ; j++ ) colors[ triangles[i][j] ] += SampleTexture( texture , tWidth , tHeight , triangles[i].data[j] , false ) , count[ triangles[i][j] ]++;
		for( int i=0 ; i<vertices.size() ; i++ ) colors[i] /= (Real)count[i];
	}
}


template< class Real >
Point2D< Real > BarycentricCoordinate( const Point2D< Real > v[3] , Point2D< Real > p )
{
	// solve for (s,t) such that:
	//		p-v[0] = s * ( v[1]-v[0] ) + t * ( v[2]-v[0] )
	//		       = s * w1 + t * w2
	Point2D< Real > w1 = v[1]-v[0] , w2 = v[2]-v[0];
	SquareMatrix< Real , 2 > M;
	M(0,0) = w1[0] , M(1,0) = w2[0];
	M(0,1) = w1[1] , M(1,1) = w2[1];
	return M.inverse() * ( p-v[0] );
}
template< class Real >
Real SquaredDistanceToLineSegment( const Point2D< Real > v[2] , Point2D< Real > p )
{
	Point2D< Real > d = v[1] - v[0];
	if     ( Point2D< Real >::Dot( p-v[0] , d )<=0 ) return ( p-v[0] ).squareNorm();
	else if( Point2D< Real >::Dot( p-v[1] , d )>=0 ) return ( p-v[1] ).squareNorm();
	else
	{
		p -= v[0];
		Real dot = Point2D< Real >::Dot( p , d );
		return p.squareNorm() - dot * dot / d.squareNorm();
	}
}
template< class Real >
Real SquaredDistanceToLineSegment( Point2D< Real > v1 , Point2D< Real > v2 , Point2D< Real > p )
{
	Point2D< Real > v[] = { v1 , v2 };
	return SquaredDistanceToLineSegment( v , p );
}
template< class Real >
Real SquaredDistanceToTriangle( const Point2D< Real > v[3] , Point2D< Real > p )
{
	Point2D< Real > b = BarycentricCoordinate( v , p );
	if( b[0]>=0 && b[1]>=0 && b[0]+b[1]<=1 ) return 0;
	else return std::min< Real >( SquaredDistanceToLineSegment( v[0] , v[1] , p  ) , std::min< Real >( SquaredDistanceToLineSegment( v[1] , v[2] , p ) , SquaredDistanceToLineSegment( v[2] , v[0] , p ) ) );
}



template< class Real >
void RasterizeTriangleInterior( Point2D< Real > v0 , Point2D< Real > v1 , Point2D< Real > v2 , int tIdx , FEM::SamplePoint< Real >* samplePoints , int width , int height )
{
	// If things are done correctly, then every rasterized pixel should have the property that if we evaluate the barycentrictrically weighted average
	// of the texture coordinates of the associated triangle, we should get back the pixel coordinates.
	Point2D< double > v[] = { Point2D< double >( v0 ) , Point2D< double >( v1 ) , Point2D< double >( v2 ) };
	for( int j=0 ; j<3 ; j++ ) v[j][0] *= (width-1) , v[j][1] *= (height-1);
	// Sort the points from highest to lowest
	int map[3];
	if( v0[1]<=v1[1] && v0[1]<=v2[1] )
	{
		map[0] = 0;
		if( v1[1]<=v2[1] ) map[1] = 1 , map[2] = 2;
		else               map[1] = 2 , map[2] = 1;
	}
	else if( v1[1]<=v0[1] && v1[1]<=v2[1] )
	{
		map[0] = 1;
		if( v0[1]<=v2[1] ) map[1] = 0 , map[2] = 2;
		else               map[1] = 2 , map[2] = 0;
	}
	else
	{
		map[0] = 2;
		if( v0[1]<=v1[1] ) map[1] = 0 , map[2] = 1;
		else               map[1] = 1 , map[2] = 0;
	}
	Point2D< double > w[] = { v[ map[0] ] , v[ map[1] ] , v[ map[2] ] };
	if( w[0][1]>w[1][1] || w[1][1]>w[2][1] ) fprintf( stderr , "[ERROR] Points are out of order\n" ) , exit( 0 );
	int yStart = (int)floor( w[0][1] ) , yEnd = (int)ceil( w[2][1] );
	yStart = std::max< int >( 0 , std::min< int >( height-1 , yStart ) );
	yEnd   = std::max< int >( 0 , std::min< int >( height-1 , yEnd   ) );

	Point2D< double > source , slopes[2];
	source = w[0] , slopes[0] = w[1]-w[0] , slopes[1] = w[2]-w[0];
	for( int y=yStart ; y<=yEnd ; y++ )
	{
		if( y>=w[1][1] ) source = w[2] , slopes[0] = w[1]-w[2] , slopes[1] = w[0]-w[2];
		if( slopes[0][1]==0 || slopes[1][1]==0 ) continue;
		// source[1] + t * slopes[i][1] = y
		// => t = ( y - source[1] ) / slopes[i][1]
		// => x = sources[0] + ( y - source[1] ) * slopes[i][0] / slopes[i][1]
		double xIntercepts[] = { source[0] + ( (double)y-source[1] ) * slopes[0][0] / slopes[0][1] , source[0] + ( (double)y-source[1] ) * slopes[1][0] / slopes[1][1] };
		int xStart , xEnd;
		if( xIntercepts[0]<=xIntercepts[1] ) xStart = (int)floor( xIntercepts[0] ) , xEnd = (int)ceil( xIntercepts[1] );
		else                                 xStart = (int)floor( xIntercepts[1] ) , xEnd = (int)ceil( xIntercepts[0] );
		xStart = std::max< int >( 0 , std::min< int >( width-1 , xStart ) );
		xEnd   = std::max< int >( 0 , std::min< int >( width-1 , xEnd   ) );
		for( int x=xStart ; x<=xEnd ; x++ )
		{
			Point2D< double > b = BarycentricCoordinate( v , Point2D< double >( (double)x , (double)y ) );
			FEM::SamplePoint< Real >& p = samplePoints[ y*width + x ];
			if( b[0]>=0 && b[1]>=0 && b[0]+b[1]<=1 )
			{
				if( p.tIdx!=-1 && b[0]>1e-12 && b[1]>1e-12 && 1-b[0]-b[1]>1e-12 ) fprintf( stderr, "[WARNING] over-writing pixel ownership (%g %g %g)\n" , b[0] , b[1] , 1-b[0]-b[1] );
				p.tIdx = tIdx , p.p = Point2D< Real >( b );
			}
		}
	}
}

template< class Real >
void RemapSamplePoint( const FEM::RiemannianMesh< Real >* mesh , ConstPointer( FEM::CoordinateXForm< Real > ) xForms , FEM::SamplePoint< Real >& p )
{
	if( p.p[0]>=0 && p.p[1]>=0 && p.p[0]+p.p[1]<=1 ) return;
	else
	{
		FEM::HermiteSamplePoint< Real > _p;
		_p.tIdx = p.tIdx , _p.p = Point2D< Real >( (Real)1./3 , (Real)1./3 ) , _p.v = p.p - _p.p;
		mesh->exp( xForms , _p );
		p = _p;
	}
}

template< class Real >
FEM::SamplePoint< Real >* GetTextureSource
(
	const FEM::RiemannianMesh< Real >* mesh ,
	ConstPointer( FEM::CoordinateXForm< Real > ) xForms ,
	const std::vector< Point2D< Real > >& triangleTextures ,
	int tWidth , int tHeight ,
	int padRadius
)
{
	FEM::SamplePoint< Real >* samplePoints = new FEM::SamplePoint< Real >[ tWidth * tHeight ];
	for( int i=0 ; i<tWidth*tHeight ; i++ ) samplePoints[i].tIdx = -1;

	for( int i=0 ; i<mesh->tCount() ; i++ ) RasterizeTriangleInterior( triangleTextures[3*i] , triangleTextures[3*i+1] , triangleTextures[3*i+2] , i , samplePoints , tWidth , tHeight );

	if( padRadius>0 )
	{
		for( int i=0 ; i<tWidth ; i++ ) for( int j=0 ; j<tHeight ; j++ )
		{
			int idx = j*tWidth + i;
			if( samplePoints[idx].tIdx==-1 )
			{
				Real d2 = -1;
				Point2D< Real > p( (Real)i/(tWidth-1) , (Real)j/(tHeight-1) );
				for( int ii=i-padRadius ; ii<=i+padRadius ; ii++ ) for( int jj=j-padRadius ; jj<=j+padRadius ; jj++ ) if( ii>=0 && ii<tWidth && jj>=0 && jj<tHeight )
				{
					int _idx = jj*tWidth + ii;
					int t = samplePoints[_idx].tIdx;
					if( t>=0 )
					{
						Point2D< Real > v[] = { triangleTextures[3*t+0] , triangleTextures[3*t+1] , triangleTextures[3*t+2] };
						Real _d2 = SquaredDistanceToTriangle( v , p );
						if( d2==-1 || _d2<d2 )
						{
							d2 = _d2;
							samplePoints[idx].tIdx = -(t+2);
							samplePoints[idx].p = BarycentricCoordinate( v , p );
						}
					}
				}
			}
		}
		for( int i=0 ; i<tWidth*tHeight ; i++ ) if( samplePoints[i].tIdx<-1 ) samplePoints[i].tIdx = -(samplePoints[i].tIdx+2);
	}

	for( int i=0 ; i<tWidth*tHeight ; i++ ) if( samplePoints[i].tIdx!=-1 ) RemapSamplePoint( mesh , xForms , samplePoints[i] );
	return samplePoints;
}


#endif // TEXTURE_INCLUDED