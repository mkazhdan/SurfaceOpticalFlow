#ifndef SURFACE_OPTICAL_FLOW_INCLUDED
#define SURFACE_OPTICAL_FLOW_INCLUDED


/////////////////
// Declaration //
/////////////////


//////////////////////
// Output functions //
template< class Real >
void OutputImage( const char* fileName , const Point3D< Real >* pixels , int width , int height , bool flipY );
//////////////////////
template< class Real >
void OutputColorMesh( const char* fileName , const std::vector< Point3D< Real > >& vertices , const std::vector< Point3D< Real > >& colors , const std::vector< TriangleIndex >& triangles , int file_type );
//////////////////////
template< class Real , int Channels >
void OutputColorMesh( const char* fileName , const std::vector< Point3D< Real > >& vertices , const std::vector< Point< Real , Channels > >& colors , const std::vector< TriangleIndex >& triangles , int file_type );
//////////////////////
template< class Real , class VectorField >
void OutputFlowMesh( const char* fileName , const std::vector< Point3D< Real > >& vertices , const VectorField& vf , const std::vector< TriangleIndex >& triangles , int file_type );
// Output functions //
//////////////////////

/////////////////////////////
// Inner-product Functions //
template< class Real >
Real SquareDifference( ConstPointer( Real ) values1 , ConstPointer( Real ) values2 , size_t size );
/////////////////////////////
template< class Real >
Real Dot( ConstPointer( Real ) values1 , ConstPointer( Real) values2 , size_t size );
/////////////////////////////
template< class Real >
Real SquareNorm( ConstPointer( Real ) values , size_t size );
// Inner-product Functions //
/////////////////////////////

template< class Real , unsigned int BasisType , int Channels >
struct FlowData
{
	FlowData( void ){ mesh = NULL; }
	~FlowData( void ){ if( mesh ) delete mesh ; mesh = NULL; }
	FEM::RiemannianMesh< Real >* mesh;
	typename FEM::TangentVectorFieldWrapper< Real , BasisType >* vf;
	std::vector< Point3D< Real > > vertices;
	std::vector< TriangleIndex > triangles;
	std::vector< Point< Real , Channels > > signals[2];
	EigenSolver< Real , typename SparseMatrix< Real , int >::RowIterator > *sSolver , *vfSolver;
	SparseMatrix< Real , int > sM , sMass , sStiffness;
	SparseMatrix< Real , int > vfM , vfMass , vfStiffness;
	SparseMatrix< Real , int > Q;
	Pointer( FEM::CoordinateXForm< Real > ) xForms;
	std::vector< Real > flowCoefficients;

	void init( void );
	Real getError( bool exponential , bool halfway ) const;
	void setFittingSystem( const std::vector< Point< Real , Channels > > resampled[2] , Pointer( Real ) b , bool halfWay );
#ifdef USE_MASS
	void combineFitnessMassAndStiffness( Real massWeight , Real smoothWeight );
#else // !USE_MASS
	void combineFitnessAndStiffness( Real smoothWeight );
#endif // USE_MASS
	void smoothSignal( std::vector< Point< Real , Channels > > smoothed[2] , Real smoothWeight , bool verbose );

	template< class V >
	void transportSignal( ConstPointer( V ) in , Pointer( V ) out , Real length , bool exponential ) const;
};

template< class Real >
struct InputGeometryData
{
	std::vector< Point3D< Real > > vertices[2];
	std::vector< Point3D< Real > > colors[2];

	template< unsigned int BasisType , int Channels > void transport( const FlowData< Real , BasisType , Channels >& flowData , Real alpha , std::vector< Point3D< Real > > outputColors[2] , bool exponential );
};
template< class Real >
struct InputTextureData
{
	std::vector< Point3D< Real > > vertices;
	int tWidth , tHeight;
	unsigned char* textures[2];
	FEM::SamplePoint< Real >* textureSource;
	std::vector< Point2D< Real > > triangleTextures;

	InputTextureData( void );
	template< unsigned int BasisType , int Channels > void transport( const FlowData< Real , BasisType , Channels >& flowData , int frames , Point3D< Real >** outputTextures[2] , bool exponential , bool nearest );
};


////////////////
// Definition //
////////////////

template< class Real >
void OutputImage( const char* fileName , const Point3D< Real >* pixels , int width , int height , bool flipY )
{
	unsigned char* _pixels = new unsigned char[ width * height * 3 ];

	for( int i=0 ; i<width ; i++ ) for( int j=0 ; j<height ; j++ ) for( int c=0 ; c<3 ; c++ )
		if( flipY ) _pixels[3*((height-1-j)*width+i)+c] = (unsigned char)std::max< int >( 0 , std::min< int >( 255 , (int)floor( pixels[j*width+i][c]+0.5 ) ) );
		else        _pixels[3*((         j)*width+i)+c] = (unsigned char)std::max< int >( 0 , std::min< int >( 255 , (int)floor( pixels[j*width+i][c]+0.5 ) ) );
	char* ext = GetFileExtension( fileName );
	if( !strcasecmp( ext , "png" ) ) PNGWriteColor( fileName , _pixels , width , height );
	else if( !strcasecmp( ext , "jpg" ) || !strcasecmp( ext , "jpeg" ) )JPEGWriteColor( fileName , _pixels , width , height , 100 );
	else fprintf( stderr , "[ERROR] Unrecognized image extension: %s\n" , ext ) , exit( 0 );
	delete[] _pixels;
}

template< class Real >
void OutputColorMesh( const char* fileName , const std::vector< Point3D< Real > >& vertices , const std::vector< Point3D< Real > >& colors , const std::vector< TriangleIndex >& triangles , int file_type )
{
	std::vector< PlyColorVertex< float > > _vertices( vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ )
	{
		_vertices[i].point = Point3D< float >( vertices[i] ) , _vertices[i].color = Point3D< float >( colors[i] );
		for( int j=0 ; j<3 ; j++ ) _vertices[i].color[j] = std::min< float >( 255.f , std::max< float >( 0.f , _vertices[i].color[j] ) );
	}
	PlyWriteTriangles( fileName , _vertices , triangles , PlyColorVertex< float >::WriteProperties , PlyColorVertex< float >::WriteComponents , file_type );
}

template< class Real , int Channels >
void OutputColorMesh( const char* fileName , const std::vector< Point3D< Real > >& vertices , const std::vector< Point< Real , Channels > >& colors , const std::vector< TriangleIndex >& triangles , int file_type )
{
	if( Channels!=3 && Channels!=6 ){ fprintf( stderr , "[WARNING] Can only output mesh for 3 and 6 channel signals\n" ) ; return; }
	std::vector< PlyColorVertex< float > > _vertices( vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ )
	{
		Point< float , Channels > color;
		_vertices[i].point = Point3D< float >( vertices[i] ) , color = Point< float , Channels >( colors[i] );
		if     ( Channels==3 ) for( int j=0 ; j<3 ; j++ ) _vertices[i].color[j] = std::min< float >( 255.f , std::max< float >( 0.f , color[j] ) );
		else if( Channels==6 ) for( int j=0 ; j<3 ; j++ ) _vertices[i].color[j] = std::min< float >( 255.f , std::max< float >( 0.f , color[j] + color[j+3] ) );
	}
	PlyWriteTriangles( fileName , _vertices , triangles , PlyColorVertex< float >::WriteProperties , PlyColorVertex< float >::WriteComponents , file_type );
}

template< class Real , class VectorField >
void OutputFlowMesh( const char* fileName , const std::vector< Point3D< Real > >& vertices , const VectorField& vf , const std::vector< TriangleIndex >& triangles , int file_type )
{
	std::vector< PlyVertex< float > > _vertices( vertices.size() );
	std::vector< PlyVFFace< float > > _triangles( triangles.size() );
	FEM::SamplePoint< Real > p;
	p.p = Point2D< Real >( (Real)1./3 , (Real)1./3 );
	for( int i=0 ; i<triangles.size() ; i++ )
	{
		p.tIdx = i;
		_triangles[i].resize( 3 );
		for( int j=0 ; j<3 ; j++ ) _triangles[i][j] = triangles[i][j];

		Point2D< Real > v = vf( p );
		_triangles[i].v = 
			( vertices[ triangles[i][1] ] - vertices[ triangles[i][0] ] ) * (float)v[0] +
			( vertices[ triangles[i][2] ] - vertices[ triangles[i][0] ] ) * (float)v[1] ;
	}
	for( int i=0 ; i<vertices.size() ; i++ ) _vertices[i].point = Point3D< float >( vertices[i] );
	PlyWritePolygons( fileName , _vertices , _triangles , PlyVertex< float >::WriteProperties , PlyVertex< float >::WriteComponents , PlyVFFace< float >::WriteProperties , PlyVFFace< float >::WriteComponents , PLY_BINARY_NATIVE );
}

template< class Real >
Real SquareDifference( ConstPointer( Real ) values1 , ConstPointer( Real ) values2 , size_t size )
{
	Real diff2 = (Real)0;
#pragma omp parallel for reduction( + : diff2 )
	for( int i=0 ; i<size ; i++ ) diff2 += ( values1[i] - values2[i] ) * ( values1[i] - values2[i] );
	return diff2;
}

template< class Real >
Real Dot( ConstPointer( Real ) values1 , ConstPointer( Real) values2 , size_t size )
{
	Real dot = (Real)0;
#pragma omp parallel for reduction( + : dot )
	for( int i=0 ; i<size ; i++ ) dot += values1[i] * values2[i];
	return dot;
}

template< class Real >
Real SquareNorm( ConstPointer( Real ) values , size_t size ){ return Dot( values , values , size ); }

//////////////
// FlowData //
//////////////
template< class Real , unsigned int BasisType , int Channels >
void FlowData< Real , BasisType , Channels >::init( void )
{
	xForms = mesh->getCoordinateXForms();

	sMass = mesh->template massMatrix< FEM::BASIS_0_WHITNEY >( );
	sStiffness = mesh->template stiffnessMatrix< FEM::BASIS_0_WHITNEY >( );
	sM.resize( sMass.rows );
#pragma omp parallel for
	for( int i=0 ; i<sMass.rows ; i++ )
	{
		sM.SetRowSize( i , sMass.rowSizes[i] );
		for( int j=0 ; j<sMass.rowSizes[i] ; j++ ) sM[i][j] = MatrixEntry< Real , int >( sMass[i][j].N , (Real)0 );
	}
	sSolver = new EigenSolverCholeskyLLt< Real , typename SparseMatrix< Real , int >::RowIterator >( sM , true );
	vfMass = mesh->template massMatrix< BasisType >( );
	vfStiffness = mesh->template stiffnessMatrix< BasisType >( );

#ifdef USE_MASS
	{
		if( vfStiffness.rows!=vfMass.rows ){ fprintf( stderr , "[ERROR] Rows don't match\n" ) ; exit(0); }

		SparseMatrix< Real , int > _vfMass;
		_vfMass.resize( vfMass.rows );
		for( unsigned int r=0 ; r<vfStiffness.rows ; r++ )
		{
			_vfMass.SetRowSize( r , vfStiffness.rowSizes[r] );
			for( unsigned int c=0 ; c<vfStiffness.rowSizes[r] ; c++ ) _vfMass[r][c].N = vfStiffness[r][c].N , _vfMass[r][c].Value = 0;

			for( unsigned int c=0 ; c<vfMass.rowSizes[r] ; c++ )
			{
				bool found = false;
				for( unsigned int _c=0 ; _c<_vfMass.rowSizes[r] ; _c++ )
					if( vfMass[r][c].N==_vfMass[r][_c].N )
					{
						found = true;
						_vfMass[r][_c].Value = vfMass[r][c].Value;
					}
				if( !found ){ fprintf( stderr , "[ERROR] Couldn't find mass matrix entry in stiffness\n" ) ; exit(0); }
			}
		}
		vfMass = _vfMass;
	}
#endif // USE_MASS

	flowCoefficients.resize( mesh->template dimension< BasisType >() , (Real)0 );

	{
		Pointer( SquareMatrix< Real , 2 > ) newTensors = AllocPointer< SquareMatrix< Real , 2 > >( triangles.size() );
#pragma omp parallel for
		for( int i=0 ; i<triangles.size() ; i++ ) newTensors[i] = SquareMatrix< Real , 2 >::Identity();
		Q = mesh->template massMatrix< BasisType >( false , newTensors );
		FreePointer( newTensors );
	}

		// Reorder the stiffness matrix to match the order of the fitting matrix
#pragma omp parallel for
	for( int i=0 ; i<Q.rows ; i++ ) for( int j=0 ; j<Q.rowSizes[i] ; j++ )
	{
		int count = 0;
		for( int k=0 ; k<vfStiffness.rowSizes[i] ; k++ ) if( Q[i][j].N==vfStiffness[i][k].N )
		{
			{
				MatrixEntry< Real , int > temp = vfStiffness[i][j];
				vfStiffness[i][j] = vfStiffness[i][k];
				vfStiffness[i][k] = temp;
			}
#ifdef USE_MASS
			{
				MatrixEntry< Real , int > temp = vfMass[i][j];
				vfMass[i][j] = vfMass[i][k];
				vfMass[i][k] = temp;
			}
#endif // USE_MASS
			count++;
		}
		if     ( count==1 ) ;
		else if( count==0 )
		{
			{
				size_t k = vfStiffness.rowSizes[i];
				vfStiffness.ResetRowSize( i , k+1 );
				vfStiffness[i][k] = vfStiffness[i][j];
				vfStiffness[i][j] = MatrixEntry< Real , int >( Q[i][j].N , 0 );
			}
#ifdef USE_MASS
			{
				size_t k = vfMass.rowSizes[i];
				vfMass.ResetRowSize( i , k+1 );
				vfMass[i][k] = vfMass[i][j];
				vfMass[i][j] = MatrixEntry< Real , int >( Q[i][j].N , 0 );
			}
#endif // USE_MASS
		}
		else fprintf( stderr , "[ERROR] FlowData::InitializeFittingSystem: %d > 1\n" , count ) , exit(0);
	}
	vfM.resize( vfStiffness.rows );
#pragma omp parallel for
	for( int i=0 ; i<vfStiffness.rows ; i++ )
	{
		vfM.SetRowSize( i , vfStiffness.rowSizes[i] );
		for( int j=0 ; j<vfStiffness.rowSizes[i] ; j++ ) vfM[i][j] = MatrixEntry< Real , int >( vfStiffness[i][j].N , (Real)0 );
	}

	if( FEM::BasisInfo< BasisType >::Singular ) vfSolver = new EigenSolverCholeskyLDLt< Real , typename SparseMatrix< Real , int >::RowIterator >( vfM , true );
	else                                        vfSolver = new EigenSolverCholeskyLLt < Real , typename SparseMatrix< Real , int >::RowIterator >( vfM , true );

	vf = new typename  FEM::TangentVectorFieldWrapper< Real , BasisType >( mesh , ( ConstPointer( Real ) )GetPointer( flowCoefficients ) , true );
}

template< class Real , unsigned int BasisType , int Channels >
Real FlowData< Real , BasisType , Channels >::getError( bool exponential , bool halfway ) const
{
	Real error = (Real)0;
	Real stepSize = (Real)( halfway ? 0.5 : 1. );
	ConstPointer( Point< Real , Channels > ) _signals[] = { GetPointer( signals[0] ) , GetPointer( signals[1] ) };
	ConstPointer( TriangleIndex ) _triangles = GetPointer( triangles );
#pragma omp parallel for reduction ( + : error )
	for( int i=0 ; i<triangles.size() ; i++ )
	{
		FEM::SamplePoint< Real > p[] = { FEM::SamplePoint< Real >( i , Point2D< Real >( (Real)1./3 , (Real)1./3 ) ) , FEM::SamplePoint< Real >( i , Point2D< Real >( (Real)1./3 , (Real)1./3 ) ) };
		FEM::SamplePoint< Real > q[] = { p[0] , p[1] };
		if( exponential ) for( int s=0 ; s<2 ; s++ )
		{
			FEM::HermiteSamplePoint< Real > _p( p[s] , (*vf)( p[s] ) * ( s==0 ? -stepSize : stepSize ) );
			mesh->exp( xForms , _p );
			p[s] = _p;
		}
		else for( int s=0 ; s<2 ; s++ ) mesh->flow( xForms , *vf , s==0 ? -stepSize : stepSize , p[s] , (Real)FLOW_STEP_SIZE );
		Real a = mesh->area( i );
		if( halfway ) error += Point< Real , Channels >::SquareNorm( mesh->template evaluateScalarField< FEM::BASIS_0_WHITNEY >( _signals[0] , p[0] ) - mesh->template evaluateScalarField< FEM::BASIS_0_WHITNEY >( _signals[1] , p[1] ) ) * a;
		else
		{
			error += Point< Real , Channels >::SquareNorm( mesh->template evaluateScalarField< FEM::BASIS_0_WHITNEY >( _signals[0] , q[0] ) - mesh->template evaluateScalarField< FEM::BASIS_0_WHITNEY >( _signals[1] , p[1] ) ) * a;
			error += Point< Real , Channels >::SquareNorm( mesh->template evaluateScalarField< FEM::BASIS_0_WHITNEY >( _signals[0] , p[0] ) - mesh->template evaluateScalarField< FEM::BASIS_0_WHITNEY >( _signals[1] , q[1] ) ) * a;
		}
	}
	return halfway ?  error : error / (Real)2.;
}

template< class Real , unsigned int BasisType , int Channels >
void FlowData< Real , BasisType , Channels >::smoothSignal( std::vector< Point< Real , Channels > > smoothed[2] , Real smoothWeight , bool verbose )
{
	double setTime;
	for( int s=0 ; s<2 ; s++ ) smoothed[s].resize( vertices.size() );
#pragma omp parallel for
	for( int i=0 ; i<sM.rows ; i++ ) for( int j=0 ; j<sM.rowSizes[i] ; j++ ) sM[i][j].Value = sMass[i][j].Value + sStiffness[i][j].Value * smoothWeight;
	{
		Timer t;
		sSolver->update( sM );
		setTime = t.elapsed();
	}
	Timer t;
	Pointer( Real ) x = AllocPointer< Real >( vertices.size() );
	Pointer( Real ) b = AllocPointer< Real >( vertices.size() );
	for( int s=0 ; s<2 ; s++ ) for( int c=0 ; c<Channels ; c++ )
	{
#pragma omp parallel for
		for( int i=0 ; i<vertices.size() ; i++ ) x[i] = signals[s][i][c];
		sMass.Multiply( x , b );
		sSolver->solve( ( ConstPointer( Real ) )b , x );
#pragma omp parallel for
		for( int i=0 ; i<vertices.size() ; i++ ) smoothed[s][i][c] = x[i];
	}
	if( verbose ) printf( "\tSet / Solved system: %.2f / %.2f (s) [%d x %d: %f]\n" , setTime , t.elapsed() , (int)sM.rows , (int)sM.rows , (Real)sM.Entries()/sM.rows );
	FreePointer( x );
	FreePointer( b );
}

template< class Real , unsigned int BasisType , int Channels >
void FlowData< Real , BasisType , Channels >::setFittingSystem( const std::vector< Point< Real , Channels > > resampled[2] , Pointer( Real ) b , bool halfWay )
{
	Q *= 0;
	memset( b , 0 , sizeof(Real) * vfM.rows );
	auto OuterProduct = []( Point2D< Real > v )
	{
		SquareMatrix< Real , 2 > m;
		for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) m(i,j) = v[i] * v[j];
		return m;
	};
	Pointer( SquareMatrix< Real , 2 > ) newTensors = AllocPointer< SquareMatrix< Real , 2 > >( triangles.size() );
#pragma omp parallel for
	for( int t=0 ; t<triangles.size() ; t++ )
	{
		SquareMatrix< Real , 2 > newTensor;
		FEM::CotangentVector< Real > cGrads[3];
		for( int c=0 ; c<Channels ; c++ )
		{
			FEM::CotangentVector< Real > grads[2];
			for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<3 ; j++ ) grads[i] += FEM::RightTriangle< Real >::CornerGradients[j] * resampled[i][ triangles[t][j] ][c];

			auto Process = [&]( FEM::CotangentVector< Real > grad )
			{
				newTensor += OuterProduct( grad );
				for( int j=0 ; j<3 ; j++ ) cGrads[j] += grad * ( resampled[0][ triangles[t][j] ][c] - resampled[1][ triangles[t][j] ][c] );
			};
			if( halfWay ) Process( ( grads[0] + grads[1] ) / 2 );
			else Process( grads[0] ) , Process( grads[1] );
		}
		Point< Real , FEM::BasisInfo< BasisType >::Coefficients > iDual = FEM::RightTriangle< Real >::template IntegrationDual< BasisType >( mesh->g(t) , ( ConstPointer( FEM::CotangentVector< Real > ) )GetPointer( cGrads , 3 ) );
		for( int j=0 ; j<FEM::BasisInfo< BasisType >::Coefficients ; j++ )
		{
			bool aligned;
			int e = mesh->template index< BasisType >( t , j , aligned );
			if( aligned )
#pragma omp atomic
				b[e] += iDual[j];
			else
#pragma omp atomic
				b[e] -= iDual[j];
		}
		newTensors[t] = newTensor;
	}

	Q = mesh->template massMatrix< BasisType >( false , newTensors );
	FreePointer( newTensors );
}
template< class Real , unsigned int BasisType , int Channels >
#ifdef USE_MASS
void FlowData< Real , BasisType , Channels >::combineFitnessMassAndStiffness( Real massWeight , Real smoothWeight )
#else // !USE_MASS
void FlowData< Real , BasisType , Channels >::combineFitnessAndStiffness( Real smoothWeight )
#endif // USE_MASS
{
#pragma omp parallel for
#ifdef USE_MASS
	for( int i=0 ; i<vfStiffness.rows ; i++ ) for( int j=0 ; j<vfStiffness.rowSizes[i] ; j++ ) vfM[i][j].Value = vfStiffness[i][j].Value * smoothWeight + vfMass[i][j].Value * massWeight;
#else // !USE_MASS
	for( int i=0 ; i<vfStiffness.rows ; i++ ) for( int j=0 ; j<vfStiffness.rowSizes[i] ; j++ ) vfM[i][j].Value = vfStiffness[i][j].Value * smoothWeight;
#endif // USE_MASS
#pragma omp parallel for
	for( int i=0 ; i<Q.rows ; i++ ) for( int j=0 ; j<Q.rowSizes[i] ; j++ ) vfM[i][j].Value += Q[i][j].Value;
}

template< class Real , unsigned int BasisType , int Channels >
template< class V >
void FlowData< Real , BasisType , Channels >::transportSignal( ConstPointer( V ) in , Pointer( V ) out , Real length , bool exponential ) const
{
	std::vector< V > _out( mesh->tCount()*3 );
	std::vector< int > counts( mesh->vCount() , 0 );

#pragma omp parallel for
	for( int i=0 ; i<mesh->vCount() ; i++ ) out[i] *= (Real)0;
#pragma omp parallel for
	for( int i=0 ; i<_out.size() ; i++ ) _out[i] *= (Real)0;

#pragma omp parallel for
	for( int i=0 ; i<mesh->tCount() ; i++ )
	{
		FEM::SamplePoint< Real > p( i , Point2D< Real >( (Real)1./3 , (Real)1./3 ) );
		if( exponential )
		{
			FEM::HermiteSamplePoint< Real > _p( p , (*vf)(p)*length );
			mesh->exp( xForms , _p );
			p = _p;
		}
		else mesh->flow( xForms , *vf , length , p , (Real)FLOW_STEP_SIZE );
		V c = mesh->template evaluateScalarField< FEM::BASIS_0_WHITNEY >( in , p );
		for( int j=0 ; j<3 ; j++ ) _out[ 3*i+j ] += c;
	}

	for( int i=0 ; i<mesh->tCount() ; i++ ) for( int j=0 ; j<3 ; j++ ) out[ mesh->triangles(i)[j] ] += _out[3*i+j] , counts[ mesh->triangles(i)[j] ]++;

#pragma omp parallel for
	for( int i=0 ; i<mesh->vCount() ; i++ ) out[i] /= (Real)counts[i];
}

///////////////////////
// InputGeometryData //
///////////////////////

template< class Real >
template< unsigned int BasisType , int Channels >
void InputGeometryData< Real >::transport( const FlowData< Real , BasisType , Channels >& flowData , Real alpha , std::vector< Point3D< Real > > outputColors[2] , bool exponential )
{
	for( int s=0 ; s<2 ; s++ )
	{
		Real length = (Real)( s==0 ? -alpha : 1.-alpha );
		outputColors[s].resize( colors[s].size() );
		flowData.transportSignal( ( ConstPointer( Point3D< Real > ) )GetPointer( colors[s] ) , GetPointer( outputColors[s] ) , length , exponential );
	}
}

template< class Real >
InputTextureData< Real >::InputTextureData( void ){ tWidth = tHeight = 0 , textures[0] = textures[1] = NULL; }

//////////////////////
// InputTextureData //
//////////////////////

template< class Real >
template< unsigned int BasisType , int Channels >
void InputTextureData< Real >::transport( const FlowData< Real , BasisType , Channels >& flowData , int frames , Point3D< Real >** outputTextures[2] , bool exponential , bool nearest )
{
	const typename FEM::TangentVectorFieldWrapper< Real , BasisType >& vf = *flowData.vf;
	for( int s=0 ; s<2 ; s++ )
	{
		Real alpha = (Real)1./(frames-1);
		Real length = (Real)( s==0 ? -alpha : alpha );
#pragma omp parallel for
		for( int i=0 ; i<tWidth*tHeight ; i++ ) if( textureSource[i].tIdx!=-1 )
		{
			FEM::SamplePoint< Real > p = textureSource[i];
			for( int f=0 ; f<frames ; f++ )
			{
				Point2D< Real > q = triangleTextures[ p.tIdx*3+0 ] * ( (Real)1. - p.p[0] - p.p[1] ) + triangleTextures[ p.tIdx*3+1 ] * p.p[0] + triangleTextures[ p.tIdx*3+2 ] * p.p[1];
				outputTextures[s][f][i] = SampleTexture( textures[s] , tWidth , tHeight , q , !nearest );
				if( exponential )
				{
					FEM::HermiteSamplePoint< Real > _p( p , vf(p)*length );
					flowData.mesh->exp( flowData.xForms , _p );
					p = _p;
				}
				else flowData.mesh->flow( flowData.xForms , vf , length , p , (Real)( FLOW_STEP_SIZE * frames ) );
			}
		}
	}
}
#endif // SURFACE_OPTICAL_FLOW_INCLUDED