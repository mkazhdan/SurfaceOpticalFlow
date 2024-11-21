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


#undef ARRAY_DEBUG
#define FOR_RELEASE

#define FLOW_STEP_SIZE 1e-3

#include "Src/PreProcessor.h"

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <algorithm>
#include <vector>
#include <Misha/cmdLineParser.h>
#include <Misha/Geometry.h>
#include <Misha/Ply.h>
#include <Misha/LinearSolvers.h>
#include <Misha/Timer.h>
#include <Misha/PNG.h>
#include <Misha/JPEG.h>
#include <Misha/FEM.h>
#include "Src/Texture.inl"
#include "Src/Subdivide.inl"
#include "SurfaceOpticalFlow.inl"

const char defaultExtension[] = "png";

cmdLineParameter< char* >
	BaseMesh( "mesh" ) ,
	FlowField( "field" ) ,
	Out( "out" ) ,
	Advected( "advected" ) ,
	Extension( "ext" );

cmdLineParameterArray< char* , 2 >
	In( "in" );
cmdLineParameter< int >
	Levels( "levels" , 7 ) ,
	SubLevels( "subLevels" , 1 ) ,
	SmoothIters( "iters" , 1 ) ,
	KeyFrames( "keyFrames" , 2 ) ,
	PadRadius( "pad" , 3 ) ,
	KrylovDimension( "search" , 1 );

cmdLineParameter< float >
	ScalarSmoothWeight( "sSmooth" , 3e-3f ) ,
	VectorFieldMassWeight( "vfMass" , 0.f ) ,
	VectorFieldSmoothWeight( "vfSmooth" , 1e6f ) ,
	SubdivideEdgeLength( "eLength" , 0.008f ) ,
	DoGWeight( "dogWeight" , 1.f ) ,
	DoGSmooth( "dogSmooth" , (float)1e-4 );

cmdLineParameter< float >
	WeightMultiplier( "wMultiply" , 0.25f );

cmdLineReadable
	Verbose( "verbose" ) ,
	ShowError( "error" ) ,
	Nearest( "nearest" ) ,
	Float( "float" ) ,
	Exponential( "exp" ) ,
	AverageError( "average" ) ,
	ShowStats( "stats" ) ,
	Whitney( "whitney" );

cmdLineReadable* params[] =
{
	&BaseMesh ,
	&In ,
	&Out ,
	&Advected ,
	&Levels ,
	&SubLevels ,
	&ScalarSmoothWeight  ,
	&VectorFieldMassWeight ,
	&VectorFieldSmoothWeight ,
	&SmoothIters ,
	&Verbose ,
	&ShowError ,
	&WeightMultiplier ,
	&FlowField ,
	&DoGWeight ,
	&DoGSmooth ,
	&KeyFrames ,
	&SubdivideEdgeLength ,
	&Nearest ,
	&PadRadius ,
	&Float ,
	&KrylovDimension ,
	&Exponential ,
	&AverageError ,
	&ShowStats ,
	&Whitney ,
	&Extension ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input source / target>\n" , In.name );
	printf( "\t[--%s <input geometry>]\n" , BaseMesh.name );
	printf( "\t[--%s <hierarchy levels>=%d]\n" , Levels.name , Levels.value );
	printf( "\t[--%s <scalar smoothing weight>=%f]\n" , ScalarSmoothWeight.name , ScalarSmoothWeight.value );
	printf( "\t[--%s <vector field mass weight>=%f]\n" , VectorFieldMassWeight.name , VectorFieldMassWeight.value );
	printf( "\t[--%s <vector field smoothing weight>=%f]\n" , VectorFieldSmoothWeight.name , VectorFieldSmoothWeight.value );
	printf( "\t[--%s <number of key-frames to output>=%d]\n" , KeyFrames.name , KeyFrames.value );
	printf( "\t[--%s <output headers>]\n" , Out.name );
	printf( "\t[--%s <output image extension>]\n" , Extension.name );
#ifndef FOR_RELEASE
	printf( "\t[--%s <advected header>]\n" , Advected.name );
	printf( "\t[--%s <flow-field>]\n" , FlowField.name );
	printf( "\t[--%s <iterations within each level>=%d]\n" , SubLevels.name , SubLevels.value );
#endif // FOR_RELEASE
	printf( "\t[--%s <padding radius>=%d]\n" , PadRadius.name , PadRadius.value );
	printf( "\t[--%s <target edge length (as a fraction of the total circumference)>=%f]\n" , SubdivideEdgeLength.name , SubdivideEdgeLength.value );
#ifndef FOR_RELEASE
	printf( "\t[--%s <weight multiplication factor>=%f]\n" , WeightMultiplier.name , WeightMultiplier.value );
#endif // FOR_RELEASE
	printf( "\t[--%s <difference of Gaussians blending weight>=%f]\n" , DoGWeight.name , DoGWeight.value );
	printf( "\t[--%s <difference of Gaussians smoothing weight>=%f]\n" , DoGSmooth.name , DoGSmooth.value );
#ifndef FOR_RELEASE
	printf( "\t[--%s <Krylov dimension>=%d]\n" , KrylovDimension.name , KrylovDimension.value );
#endif // FOR_RELEASE
	printf( "\t[--%s]\n" , Whitney.name );
#ifndef FOR_RELEASE
	printf( "\t[--%s]\n" , Float.name );
	printf( "\t[--%s]\n" , Exponential.name );
	printf( "\t[--%s]\n" , AverageError.name );
	printf( "\t[--%s]\n" , Nearest.name );
	printf( "\t[--%s]\n" , ShowStats.name );
#endif // FOR_RELEASE
	printf( "\t[--%s]\n" , ShowError.name );
	printf( "\t[--%s]\n" , Verbose.name );
}

template< class Real >
struct PointData
{
	int count;
	Point3D< Real > p;
};

template< class Real >
unsigned long long PointKey( Point3D< Real > p )
{
	for( int j=0 ; j<3 ; j++ ) p[j] = std::max< Real >( 0 , std::min< Real >( 1 , p[j] ) );
	unsigned long long x , y , z;
	x = (unsigned long long)( p[0] * ( 1<<21 ) );
	y = (unsigned long long)( p[1] * ( 1<<21 ) );
	z = (unsigned long long)( p[2] * ( 1<<21 ) );
	return ( x ) | ( y<<21 ) | ( z<<42 );
}

bool IsImage( const char* fileName )
{
	char* ext = GetFileExtension( In.values[0] );
	if( !strcasecmp( ext , "png" ) || !strcasecmp( ext , "jpg" ) || !strcasecmp( ext , "jpeg" ) )
	{
		delete[] ext;
		return true;
	}
	else
	{
		delete[] ext;
		return false;
	}
}

template< class Real > struct CorrectionStats{ Real startError[2] , endError[2] , mass , stiffness; };

// Given a current estimation of the low-frequency flow, estimate the high-frequency correction flow.
template< unsigned int BasisType , unsigned int SearchDim , class Real , int Channels >
CorrectionStats< Real > EstimateCorrectionFlow
(
	FlowData< Real , BasisType , Channels >& flowData ,
	Real scalarSmoothWeight ,
	Real vectorMassWeight ,
	Real vectorSmoothWeight
)
{
	CorrectionStats< Real > stats;
	std::vector< Point< Real , Channels > > smoothed[2] , resampled[2];

	const std::vector< Point3D< Real > >& vertices = flowData.vertices;
	const std::vector< TriangleIndex >& triangles = flowData.triangles;
	ConstPointer( Point3D< Real > ) _vertices = ( ConstPointer( Point3D< Real > ) )GetPointer( vertices );

	// Smooth the signal and resample to the midway point
	{
		flowData.smoothSignal( smoothed , scalarSmoothWeight , Verbose.set );
		for( int s=0 ; s<2 ; s++ )
		{
			resampled[s].resize( smoothed[s].size() );
			flowData.transportSignal( ( ConstPointer( Point< Real , Channels > ) )GetPointer( smoothed[s] ) , GetPointer( resampled[s] ) , (Real)( s==0 ? -0.5 : 0.5 ) , Exponential.set );
		}
	}

	// Solve for the flow field
	{
		Timer t;

		Pointer( Real ) x[SearchDim];
		for( int i=0 ; i<SearchDim ; i++ ) x[i] = AllocPointer< Real >( flowData.vfM.rows );
		Pointer( Real ) _b = AllocPointer< Real >( flowData.vfM.rows );
		Pointer( Real ) b = AllocPointer< Real >( flowData.vfM.rows );

		// Compute the linear system and constraints for the linearized fitting
		flowData.setFittingSystem( resampled , b , !AverageError.set );

		// Rescale relative to the stiffness matrix (whose norm should be consistent for hat and whitney bases)
		{
			Real scale = (Real)sqrt( flowData.vfStiffness.SquareNorm() / flowData.Q.SquareNorm() );
			flowData.Q *= scale;
			for( int i=0 ; i<flowData.Q.rows ; i++ ) b[i] *= scale;
		}

		// Combine the fitting term with the smoothness term
		flowData.combineFitnessMassAndStiffness( vectorMassWeight , vectorSmoothWeight );

		double setTime = t.elapsed();

		// Solve the linear system(s)
		{
			double updateTime , solveTime;
			// Update the numerical part of the Cholesky factorixation
			{
				Timer t;
				flowData.vfSolver->update( flowData.vfM );
				updateTime = t.elapsed();
			}
			// Solve the linear system(s)
			{
				Timer t;
				memcpy( _b , b , sizeof(Real) * flowData.Q.rows );
				for( int i=0 ; i<SearchDim ; i++ )
				{
					flowData.vfSolver->solve( ( ConstPointer( Real ) )_b , x[i] );
					if( i<SearchDim-1 ) flowData.vfMass.Multiply( x[i] , _b );
				}
				solveTime = t.elapsed();
			}
			if( Verbose.set ) printf( "\tSet / Updated / Solved system: %.2f / %.2f / %.2f (s) [%d x %d: %f]\n" , setTime , updateTime , solveTime , (int)flowData.vfM.rows , (int)flowData.vfM.rows , (Real)flowData.vfM.Entries()/flowData.vfM.rows );
		}

		if( ShowError.set || ShowStats.set ) stats.startError[0] = flowData.getError( Exponential.set , false ) , stats.startError[1] = flowData.getError( Exponential.set , true );

		// Solve for the linear combination of the Krylov vectors that best solves the fitting problem (without the smoothness constraint)
		{
			SquareMatrix< Real , SearchDim > M;
			Point< Real , SearchDim > V;
			for( int i=0 ; i<SearchDim ; i++ )
			{
				{
					Real v = 0;
#pragma omp parallel for reduction( + : v )
					for( int k=0 ; k<flowData.Q.rows ; k++ ) v += x[i][k] * b[k];
					V[i] = v;
				}
				for( int j=0 ; j<=i ; j++ )
				{
					Real m = 0;
#pragma omp parallel for reduction( + : m )
					for( int k=0 ; k<flowData.Q.rows ; k++ ) for( int l=0 ; l<flowData.Q.rowSizes[k] ; l++ ) m += x[i][k] * x[j][ flowData.Q[k][l].N ] * flowData.Q[k][l].Value;
					M( i , j ) = M( j , i ) = m;
				}
			}
			V = M.inverse() * V;
#pragma omp parallel for
			for( int i=0 ; i<flowData.Q.rows ; i++ )
			{
				Real v = 0;
				for( int j=0 ; j<SearchDim ; j++ ) v += x[j][i] * V[j];
				flowData.flowCoefficients[i] += v;
			}
			if( ShowError.set || ShowStats.set ) stats.endError[0] = flowData.getError( Exponential.set , false ) , stats.endError[1] = flowData.getError( Exponential.set , true );
		}
		for( int i=0 ; i<SearchDim ; i++ ) FreePointer( x[i] );
		FreePointer( _b );
		FreePointer( b );
	}

	if( Verbose.set || ShowStats.set )
	{
		ConstPointer( Real ) x = ( ConstPointer( Real ) )GetPointer( flowData.flowCoefficients );
		Pointer( Real ) Mx = AllocPointer< Real >( flowData.vfMass.rows );
		flowData.vfMass.Multiply( x , Mx );
		stats.mass = Dot( ( ConstPointer( Real ) )Mx , ( ConstPointer( Real ) )x , flowData.vfMass.rows );
		flowData.vfStiffness.Multiply( x , Mx );
		stats.stiffness = Dot( ( ConstPointer( Real ) )Mx , ( ConstPointer( Real ) )x , flowData.vfMass.rows );
		FreePointer( Mx );
	}
	if( Verbose.set ) printf( "\t\tFlow-field mass/stiffness: %g / %g -> %g\n" , sqrt(stats.mass) , sqrt(stats.stiffness) , sqrt(stats.stiffness/stats.mass) );
	if( ShowError.set ) printf( "\t[ERROR] [%g %g] -> [%g %g]\n" , stats.startError[0] , stats.startError[1] , stats.endError[0] , stats.endError[1] );

	return stats;
}

template< unsigned int BasisType , unsigned int SearchDim , class Real , int Channels >
int __Execute( void )
{
	FlowData< Real , BasisType , Channels > flowData;
	std::vector< Point3D< Real > >& vertices = flowData.vertices;
	std::vector< TriangleIndex >& triangles = flowData.triangles;

	bool processTexture;
	InputTextureData< Real > inputTextureData;
	InputGeometryData< Real > inputGeometryData;

	int file_type;
	// Process Geometry+Textures
	if( BaseMesh.set || ( IsImage( In.values[0] ) && IsImage( In.values[1] ) ) )
	{
		processTexture = true;

		{
			int w , h;
			char* ext = GetFileExtension( In.values[0] );
			if( !strcasecmp( ext , "png" ) ) inputTextureData.textures[0] = PNGReadColor( In.values[0] , w , h );
			else if( !strcasecmp( ext , "jpg" ) || !strcasecmp( ext , "jpeg" ) ) inputTextureData.textures[0] = JPEGReadColor( In.values[0] , w , h );
			else{ fprintf( stderr , "[ERROR] Unrecognized image extension: %s\n" , ext ) ; return EXIT_FAILURE; }
			delete[] ext;
			inputTextureData.tWidth = w , inputTextureData.tHeight = h;
		}
		{
			int w , h;
			char* ext = GetFileExtension( In.values[1] );
			if( !strcasecmp( ext , "png" ) ) inputTextureData.textures[1] = PNGReadColor( In.values[1] , w , h );
			else if( !strcasecmp( ext , "jpg" ) || !strcasecmp( ext , "jpeg" ) ) inputTextureData.textures[1] = JPEGReadColor( In.values[1] , w , h );
			else{ fprintf( stderr , "[ERROR] Unrecognized image extension: %s\n" , ext ) ; return EXIT_FAILURE; }
			delete[] ext;
			if( inputTextureData.tWidth!=w || inputTextureData.tHeight!=h ){ fprintf( stderr , "[ERROR] Texture resolutions don't match: %d x %d != %d x %d\n" , inputTextureData.tWidth , inputTextureData.tHeight , w , h ) ; return EXIT_FAILURE; }
		}

		std::vector< TriangleIndexWithData< Point2D< Real > > > texturedTriangles;
		if( BaseMesh.set )
		{
			std::vector< PlyTextureVertex< float > > _vertices;
			bool hasVertexTexture , hasFaceTexture;
			{
				std::vector< PlyTexturedFace< float > > _faces;
				bool vertexFlags[ PlyTextureVertex< float >::ReadComponents ] , faceFlags[ PlyTexturedFace< float >::ReadComponents ];
				PlyReadPolygons( BaseMesh.value , _vertices , _faces , PlyTextureVertex< float >::ReadProperties , vertexFlags , PlyTextureVertex< float >::ReadComponents , PlyTexturedFace< float >::ReadProperties , faceFlags , PlyTexturedFace< float >::ReadComponents , file_type );
				hasVertexTexture = ( vertexFlags[3] && vertexFlags[4] ) || ( vertexFlags[5] && vertexFlags[6] ) , hasFaceTexture = faceFlags[1];

				MinimalAreaTriangulation< double > MAT;
				std::vector< Point3D< double > > face;
				std::vector< TriangleIndex > triangles;
				for( size_t i=0 ; i<_faces.size() ; i++ )
				{
					if( hasFaceTexture && _faces[i].nr_uv_coordinates!=2*_faces[i].size() ) fprintf( stderr , "[ERROR] Bad textured face: vertices: %d ; texture-coordinates: %d\n" , _faces[i].nr_vertices , _faces[i].nr_uv_coordinates ) , exit( 0 );

					face.resize( _faces[i].size( ) );
					MAT.GetTriangulation( face , triangles );
					for( size_t j=0 ; j<triangles.size() ; j++ )
					{
						TriangleIndexWithData< Point2D< Real > > triangle;
						for( int k=0 ; k<3 ; k++ ) triangle[k] = _faces[i][ triangles[j][k] ];
						if( hasFaceTexture ) for( int k=0 ; k<3 ; k++ ) triangle.data[k] = Point2D< Real >( _faces[i].texture( triangles[j][k] ) );
						texturedTriangles.push_back( triangle );
					}
				}
			}

			if( hasFaceTexture )
			{
				vertices.resize( _vertices.size() );
				for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = Point3D< Real >( _vertices[i].point );
			}
			else if( hasVertexTexture )
			{
				// Copy the the texture coordinates from the vertices to the corners
				for( int i=0 ; i<texturedTriangles.size() ; i++ ) for( int j=0 ; j<3 ; j++ ) texturedTriangles[i].data[j] = Point2D< Real >( _vertices[ texturedTriangles[i][j] ].texture[0] , _vertices[ texturedTriangles[i][j] ].texture[1] );

				// Get the bounding box of the vertex set
				Point3D< float > min , max;
				min = max = _vertices[0].point;
				for( int i=0 ; i<_vertices.size() ; i++ ) for( int j=0 ; j<3 ; j++ ) min[j] = std::min< float >( min[j] , _vertices[i].point[j] ) , max[j] = std::max< float >( max[j] , _vertices[i].point[j] );

				std::unordered_map< unsigned long long , int > vertexMap;
				std::vector< int > _vertexMap( _vertices.size() );
				std::vector< PointData< float > > pointData;
				for( int i=0 ; i<_vertices.size() ; i++ )
				{
					Point3D< float > p;
					for( int j=0 ; j<3 ; j++ ) p[j] = ( _vertices[i].point[j] - min[j] ) / ( max[j] - min[j] );
					unsigned long long key = PointKey( p );
					std::unordered_map< unsigned long long , int >::iterator iter = vertexMap.find( key );
					int idx;
					if( iter==vertexMap.end() )
						{
							PointData< float > pd;
						pd.count = 1;
						pd.p = _vertices[i].point;
						idx = (int)pointData.size();
						vertexMap[ key ] = idx;
						pointData.push_back( pd );
					}
					else
					{
						idx = iter->second;
						PointData< float >& pd = pointData[idx];
						pd.count++;
						pd.p += _vertices[i].point;
					}
					_vertexMap[i] = idx;
				}
				vertices.resize( pointData.size() );
				for( int i=0 ; i<pointData.size() ; i++ ) vertices[i] = Point3D< Real >( pointData[i].p / (float)pointData[i].count );
				for( int i=0 ; i<texturedTriangles.size() ; i++ ) for( int j=0 ; j<3 ; j++ ) texturedTriangles[i][j] = _vertexMap[ texturedTriangles[i][j] ];
			}
			else{ fprintf( stderr , "[ERROR] Mesh does not have either vertex or corner texture\n"  ) ; return EXIT_FAILURE; }
		}
		else
		{
			Point2D< Real > tCoordinates[] = { Point2D< Real >( 0 , 0 ) , Point2D< Real >( 1 , 0 ) , Point2D< Real >( 1 , 1 ) , Point2D< Real >( 0 , 1 ) };
			vertices.resize( 4 );
			for( int i=0 ; i<4 ; i++ ) vertices[i] = Point3D< Real >( tCoordinates[i][0] * inputTextureData.tWidth , tCoordinates[i][1] * inputTextureData.tHeight , 0 );
			texturedTriangles.resize( 2 );
			texturedTriangles[0][0] = 0 , texturedTriangles[0][1] =  1 , texturedTriangles[0][2] = 2;
			texturedTriangles[1][0] = 0 , texturedTriangles[1][1] =  2 , texturedTriangles[1][2] = 3;
			for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<3 ; j++ ) texturedTriangles[i].data[j] = tCoordinates[ texturedTriangles[i][j] ];
		}
		Real length = 0;
#pragma omp parallel for reduction( + : length )
		for( int i=0 ; i<texturedTriangles.size() ; i++ )
		{
			const TriangleIndexWithData< Point2D< Real > >& t = texturedTriangles[i];
			length += (Real)Length( Point3D< Real >::CrossProduct( vertices[ t[1] ] - vertices[ t[0] ] , vertices[ t[2] ] - vertices[ t[0] ] ) ) / 2;
		}
		length = (Real)sqrt( length );
		SubdivideEdgeLength.value *= (float)length;
		if( SubdivideEdgeLength.value>0 ) Subdivide( texturedTriangles , vertices , (Real)SubdivideEdgeLength.value );
		triangles.resize( texturedTriangles.size() );
		inputTextureData.triangleTextures.resize( texturedTriangles.size()*3 );
		for( int i=0 ; i<texturedTriangles.size() ; i++ )
		{
			triangles[i] = texturedTriangles[i];
			for( int j=0 ; j<3 ; j++ ) inputTextureData.triangleTextures[3*i+j] = texturedTriangles[i].data[j];
		}

		for( int s=0 ; s<2 ; s++ )
		{
			std::vector< Point3D< Real > > signal;
			SampleTextureToVertices( texturedTriangles , vertices , inputTextureData.textures[s] , inputTextureData.tWidth , inputTextureData.tHeight , signal , !Nearest.set );
			flowData.signals[s].resize( signal.size() );
			for( int i=0 ; i<signal.size() ; i++ ) for( int c=0 ; c<3 ; c++ ) flowData.signals[s][i][c] = signal[i][c];
		}
	}
	// Process Geometry+Colors
	else
	{
		processTexture = false;
		std::vector< TriangleIndex > _triangles[2];
		std::vector< PlyColorVertex< float > > _vertices[2];
		PlyReadTriangles( In.values[0] , _vertices[0] , _triangles[0] , PlyColorVertex< float >::ReadProperties , NULL , PlyColorVertex< float >::ReadComponents , file_type );
		PlyReadTriangles( In.values[1] , _vertices[1] , _triangles[1] , PlyColorVertex< float >::ReadProperties , NULL , PlyColorVertex< float >::ReadComponents , file_type );
		if( _vertices[0].size()!=_vertices[1].size() ) { fprintf( stderr , "[ERROR] Vertex counts differ: %d != %d\n" , (int)_vertices[0].size() , (int)_vertices[1].size() ) ; return EXIT_FAILURE; }
		if( _triangles[0].size()!=_triangles[1].size() ) { fprintf( stderr , "[ERROR] Different number of triangles in meshes: %d != %d\n" , (int)_triangles[0].size() , (int)_triangles[1].size() ) ; return EXIT_FAILURE; }
		triangles = _triangles[0];
		vertices.resize( _vertices[0].size() );
		for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = Point3D< Real >( _vertices[0][i] ) * ( (Real)0.5 ) + Point3D< Real >( _vertices[1][i] ) * ( (Real)0.5 );

		for( int i=0 ; i<triangles.size() ; i++ ) for( int j=0 ; j<3 ; j++ ) if( _triangles[0][i][j]!=_triangles[1][i][j] )
		{
			fprintf( stderr , "[ERROR] Triangle indices don't match: [%d,%d] %d != %d\n" , i , j , _triangles[0][i][j] , _triangles[1][i][j] );
			return EXIT_FAILURE;
		}
		for( int s=0 ; s<2 ; s++ )
		{
			inputGeometryData.colors[s].resize( _vertices[s].size() );
			inputGeometryData.vertices[s].resize( _vertices[s].size() );
			flowData.signals[s].resize( _vertices[s].size() );
			for( int i=0 ; i<_vertices[s].size() ; i++ ) inputGeometryData.colors[s][i] = Point3D< Real >( _vertices[s][i].color ) , inputGeometryData.vertices[s][i] = Point3D< Real >( _vertices[s][i].point );
			for( int i=0 ; i<_vertices[s].size() ; i++ ) for( int c=0 ; c<3 ; c++ ) flowData.signals[s][i][c] = inputGeometryData.colors[s][i][c];
		}
	}
	if( Verbose.set ) printf( "Vertices / Triangles: %d / %d\n" , (int)vertices.size() , (int)triangles.size() );

	ConstPointer( Point3D< Real > ) _vertices = ( ConstPointer( Point3D< Real > ) )GetPointer( vertices );

	// Set the mesh
	{
		flowData.mesh = new FEM::RiemannianMesh< Real >( GetPointer( triangles ) , triangles.size() );
		flowData.mesh->setMetricFromEmbedding( _vertices );
		flowData.mesh->makeUnitArea();
	}
	flowData.init();
	if( processTexture ) inputTextureData.textureSource = GetTextureSource( flowData.mesh , ( ConstPointer( FEM::CoordinateXForm< Real > ) )flowData.xForms , inputTextureData.triangleTextures , inputTextureData.tWidth , inputTextureData.tHeight , PadRadius.value );

	// Set the comparison signals
	if( DoGWeight.value>0 )
	{
		Timer t;
		Real weight = (Real)DoGSmooth.value;
		Pointer( Real ) x = AllocPointer< Real >( vertices.size() );
		Pointer( Real ) b = AllocPointer< Real >( vertices.size() );
#pragma omp parallel for
		for( int i=0 ; i<flowData.sM.rows ; i++ ) for( int j=0 ; j<flowData.sM.rowSizes[i] ; j++ ) flowData.sM[i][j].Value = flowData.sMass[i][j].Value + flowData.sStiffness[i][j].Value * weight;
		flowData.sSolver->update( flowData.sM );

		// Smooth and subtract off
		for( int s=0 ; s<2 ; s++ ) for( int c=0 ; c<3 ; c++ )
		{
			// \int [ f(x) - \int f(x) ]^2 = \int f^2(x) + [ \int f(x) ]^2 - 2 [ \int f(x) ]^2 = \int f^2(x) - [ \int f(x) ]^2

#pragma omp parallel for
			for( int i=0 ; i<vertices.size() ; i++ ) x[i] = flowData.signals[s][i][c];
			flowData.sMass.Multiply( x , b );
			Real oldAvg = flowData.mesh->getIntegral( x );
			Real oldVar = Dot( ( ConstPointer( Real ) )x , ( ConstPointer( Real ) )b , vertices.size() ) - oldAvg * oldAvg;

			flowData.sSolver->solve( ( ConstPointer( Real ) )b , x );

#pragma omp parallel for
			for( int i=0 ; i<vertices.size() ; i++ ) x[i] = flowData.signals[s][i][c] - x[i];
			flowData.sMass.Multiply( x , b );
			Real newAvg = flowData.mesh->getIntegral( x );
			Real newVar = Dot( ( ConstPointer( Real ) )x , ( ConstPointer( Real ) )b , vertices.size() ) - newAvg * newAvg;

			Real scale = (Real)sqrt( oldVar / newVar );
			if( Channels==6 )
#pragma omp parallel for
				for( int i=0 ; i<vertices.size() ; i++ ) flowData.signals[s][i][c+3] = ( x[i] - newAvg ) * scale + oldAvg;
			else if( Channels==3 )
				for( int i=0 ; i<vertices.size() ; i++ ) flowData.signals[s][i][c] = ( x[i] - newAvg ) * scale + oldAvg;
		}
		if( Channels==6 ) for( int s=0 ; s<2 ; s++ ) for( int c=0 ; c<3 ; c++ )
#pragma omp parallel for
			for( int i=0 ; i<vertices.size() ; i++ ) flowData.signals[s][i][c] *= (Real)( 1.-DoGWeight.value ) , flowData.signals[s][i][c+3] *= (Real)DoGWeight.value;
		if( Verbose.set ) printf( "Set comparison values: %.2f (s)\n" , t.elapsed() );
	}
	
	{
		Real scalarSmoothWeight = (Real)ScalarSmoothWeight.value;
		Real vectorFieldMassWeight = (Real)VectorFieldMassWeight.value;
		Real vectorFieldSmoothWeight = (Real)VectorFieldSmoothWeight.value;
		{
			Timer t;
			CorrectionStats< Real > stats;
			// Iterate over the hierarchy levels (coarse-to-fine)
			for( int i=0 ; i<Levels.value ; i++ , scalarSmoothWeight *= WeightMultiplier.value , vectorFieldMassWeight *= WeightMultiplier.value , vectorFieldSmoothWeight *= WeightMultiplier.value )
			{
				Timer t;
				// Perform the correction(s) within the level
				{
					CorrectionStats< Real > _stats;
					for( int j=0 ; j<SubLevels.value ; j++ ) _stats = EstimateCorrectionFlow< BasisType , SearchDim , Real ,  Channels >( flowData , scalarSmoothWeight , vectorFieldMassWeight , vectorFieldSmoothWeight );
					if( i==0 ) stats = _stats;
					else stats.endError[0] = _stats.endError[0] , stats.endError[1] = _stats.endError[1] , stats.mass = _stats.mass , stats.stiffness = _stats.stiffness;
				}
				if( Verbose.set ) printf( "Got flow[%d]: %.2f (s)\n" , i , t.elapsed() );
			}

			if( Verbose.set ) printf( "Got flow: %.2f (s)\n" , t.elapsed() );

			if( ShowStats.set )
			{
				printf( "E: %f -> %f , %f\n" , stats.startError[0] , stats.endError[0] , stats.endError[1] );
				printf( "S: %f\n" , sqrt( stats.stiffness / stats.mass ) );
			}
		}

		if( FlowField.set ) OutputFlowMesh( FlowField.value , vertices, *flowData.vf , triangles , file_type );

		// Generate the advected/interpolated signals
		if( Advected.set || Out.set )
		{
			Timer t;
			char fileNames[2][512] , fileName[512];
			if( processTexture )
			{
				Point3D< Real > **_textures[2] , **_texture;
				for( int s=0 ; s<2 ; s++ )
				{
					_textures[s] = new Point3D< Real >*[KeyFrames.value];
					for( int k=0 ; k<KeyFrames.value ; k++ ) _textures[s][k] = new Point3D< Real >[ inputTextureData.tWidth * inputTextureData.tHeight ];
				}
				_texture = new Point3D< Real >*[KeyFrames.value];
				for( int k=0 ; k<KeyFrames.value ; k++ ) _texture[k] = new Point3D< Real >[ inputTextureData.tWidth * inputTextureData.tHeight ];

				inputTextureData.template transport< BasisType , Channels >( flowData , KeyFrames.value , _textures , Exponential.set , Nearest.set );
#pragma omp parallel for
				for( int k=0 ; k<KeyFrames.value ; k++ )
				{
					int _k = KeyFrames.value-1-k;
					Real alpha = (Real)k / ( KeyFrames.value-1 );
					for( int i=0 ; i<inputTextureData.tWidth*inputTextureData.tHeight ; i++ ) _texture[k][i] = _textures[0][k][i] * ( (Real)(1.-alpha) ) + _textures[1][_k][i] * ( (Real)alpha );
					char fileNames[2][512] , fileName[512];
					if( Advected.set )
					{
						sprintf( fileNames[0] , "%s.S.%d.%s" , Advected.value , k , Extension.set ? Extension.value : defaultExtension );
						sprintf( fileNames[1] , "%s.T.%d.%s" , Advected.value , k , Extension.set ? Extension.value : defaultExtension );
						OutputImage( fileNames[0] , _textures[0][ k] , inputTextureData.tWidth , inputTextureData.tHeight , true );
						OutputImage( fileNames[1] , _textures[1][_k] , inputTextureData.tWidth , inputTextureData.tHeight , true );
					}
					if( Out.set )
					{
						sprintf( fileName , "%s.%d.%s" , Out.value , k , Extension.set ? Extension.value : defaultExtension );
						OutputImage( fileName , _texture[k] , inputTextureData.tWidth , inputTextureData.tHeight , true );
					}
				}
				for( int s=0 ; s<2 ; s++ )
				{
					for( int k=0 ; k<KeyFrames.value ; k++ ) delete[] _textures[s][k];
					delete[] _textures[s];
				}
				for( int k=0 ; k<KeyFrames.value ; k++ ) delete[] _texture[k];
				delete[] _texture;
			}
			else
			{
				std::vector< Point3D< Real > > advectedColors[2] , blendedColors( vertices.size() );
				for( int k=0 ; k<KeyFrames.value ; k++ )
				{
					Real alpha = (Real)k / ( KeyFrames.value-1 );
					inputGeometryData.template transport< BasisType , Channels >( flowData , alpha , advectedColors , Exponential.set );
#pragma omp parallel for
					for( int i=0 ; i<blendedColors.size() ; i++ ) blendedColors[i] = advectedColors[0][i] * ( (Real)(1.-alpha) ) + advectedColors[1][i] * ( (Real)alpha );
					if( Advected.set )
					{
						sprintf( fileNames[0] , "%s.S.%d.ply" , Advected.value , k ) , sprintf( fileNames[1] , "%s.T.%d.ply" , Advected.value , k );
						for( int s=0 ; s<2 ; s++ ) OutputColorMesh( fileNames[s] , vertices , advectedColors[s] , triangles , file_type );
					}
					if( Out.set )
					{
						sprintf( fileName , "%s.%d.ply" , Out.value , k );
						OutputColorMesh( fileName , vertices , blendedColors , triangles , file_type );
					}
				}
			}
			if( Verbose.set ) printf( "Resampled key-frames: %.2f(s)\n" , t.elapsed() );
		}
	}
	return EXIT_SUCCESS;
}
template< unsigned int SearchDim , class Real , int Channels >
int _Execute( void )
{
	if( Whitney.set ) return __Execute< FEM::BASIS_1_WHITNEY    , SearchDim , Real , Channels >( );
	else              return __Execute< FEM::BASIS_1_CONFORMING , SearchDim , Real , Channels >( );
}

template< class Real , int Channels >
int Execute( void )
{
	if( KrylovDimension.value<1 ) fprintf( stderr , "[WARNING] Setting search dimension: %d -> 1\n" , KrylovDimension.value ) , KrylovDimension.value = 1;
	switch( KrylovDimension.value )
	{
	case  1: return _Execute<  1 , Real , Channels >();
	case  2: return _Execute<  2 , Real , Channels >();
	case  4: return _Execute<  4 , Real , Channels >();
	case  8: return _Execute<  8 , Real , Channels >();
	case 16: return _Execute< 16 , Real , Channels >();
	case 32: return _Execute< 32 , Real , Channels >();
	default: fprintf( stderr , "[ERROR] Only search dimensions 1, 2, 4, 8, 16, and 32 supported\n" ) , exit( 0 );
	}
}
int main( int argc , char* argv[] )
{
	cmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}

	Timer t;
	int ret;
	DoGWeight.value = std::min< float >( 1.f , std::max< float >( 0.f , DoGWeight.value ) );
	if( Float.set )
	{
		if( DoGWeight.value>0 && DoGWeight.value<1 ) ret = Execute< float , 6 >();
		else                                         ret = Execute< float , 3 >();
	}
	else
	{
		if( DoGWeight.value>0 && DoGWeight.value<1 ) ret = Execute< double , 6 >();
		else                                         ret = Execute< double , 3 >();
	}
	if( Verbose.set ) printf( "Total running time: %.2f (s)\n" , t.elapsed() );
	return ret;
}
