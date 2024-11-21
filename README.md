<CENTER><H2>Surface Optical Flow (Version 1.0)</H2></CENTER>
<CENTER>
<A HREF="#LINKS">links</A>
<A HREF="#DESCRIPTION">description</A>
<A HREF="#EXECUTABLES">executables</A>
<A HREF="#CHANGES">changes</A>
</CENTER>
<HR>
<A NAME="LINKS"><B>LINKS</B></A><br>
<A href="https://www.cs.jhu.edu/~misha/MyPapers/SIG16.pdf">SIGGRAPH 2016 Paper</A><br>
<A HREF="https://www.cs.jhu.edu/~misha/Code/SurfaceOpticalFlow/SurfaceOpticalFlow.x64.zip">Windows (x64) Executables</A><BR>
<A href="https://www.cs.jhu.edu/~misha/Code/SurfaceOpticalFlow/SurfaceOpticalFlow.zip">Source Code</A><br> <!--<A HREF="https://github.com/mkazhdan/DMG">GitHub Repository</A><BR>-->
<!--
(Older Versions:
<A href="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version7.0/">V7.0</A>,
<A href="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version6.13a/">V6.13a</A>,
<A href="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version6.13/">V6.13</A>,
<A href="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version6.12/">V6.12</A>,
<A href="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version6.11/">V6.11</A>,
<A href="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version6.1/">V6.1</A>,
<A href="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version6/">V6</A>,
<A href="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version5.71/">V5.71</A>,
<A href="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version5.6/">V5.6</A>,
<A href="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version5.5a/">V5.5a</A>,
<A HREF="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version5.1/">V5.1</A>,
<A HREF="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version5/">V5</A>,
<A HREF="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version4.51/">V4.51</A>,
<A HREF="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version4.5/">V4.5</A>,
<A HREF="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version4/">V4</A>,
<A HREF="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version3/">V3</A>,
<A HREF="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version2/">V2</A>,
<A HREF="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version1/">V1</A>)
-->
<br>


<HR>
<A NAME="DESCRIPTION"><B>DESCRIPTION</B></A><br>
<UL>
The code takes source and target textures over a base geometry, estimates the flow field taking the source to the target and the target to the source, and outputs the in-between frames obtained by advecting the source and target along the flow (forward for the source and backward for the target) and cross-disovling the advecte results.

<P> To compute the flow field, the code performs a hierarchical solve, to find the vector field taking the source signal to the target. In the implementation, the signal can be defined by the texture colors, by the difference of Gaussians (DOG) of the texture colors, or by a weighted combination of the two.<br>
By default, the code uses the difference of Gaussians to provide a solution that is more robust to changes in lighting between the source and target.

<P> In general a base mesh in <A HREF="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format (<B>--mesh</B>) with texture coordinates in either PNG or JPEG format (<B>--in</B>) should be provided, with texture coordinates specified either at vertices or at triangle corners. In the former case, the code transforms the texture coordinates to values at corners and merges duplicate vertices.<br>
If no base geometry is specified, the code constructs a single rectangle, with the aspect ratio of the input images, essentially performing optical flow in image-space.<br>
Alternatively, colors can be specified per-vertex, in which case the two input files (<B>--in</B>) should be in <A HREF="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format and should represent the same geometry with different color values at the vertices.
  

<P> To compute the optical flow field, the code needs to solve a set of sparse linear systems of equations. The implementation supports the use of either <A HREF="https://eigen.tuxfamily.org/">Eigen</A> or <A HREF="https://faculty.cse.tamu.edu/davis/suitesparse.html">Cholmod</A>. When compiling, the type of solver is specified by defining either the <I>USE_EIGEN</I> or <I>USE_CHOLMOD</I> flag. (This is supported by the two different Makefiles / Visual Studio Projects, and the two different windows executable correspond to implementations with the two different solver libraries.)<br>
As the solution of the linear system accounts for a lion's share of the computation, we recommend using an efficient solver such as the one backed by Intel's <A HREF="https://software.intel.com/en-us/intel-mkl/">MKL</A> library. (Assuming that the libraries exist and are in a library where the linker can find them, Eigen's MKL support can be enabled by defining the <I>EIGEN_USE_MKL_ALL</I> flag.)
</UL>

<HR>
<A NAME="EXECUTABLES">
<B>Surface Optical Flow:</B><BR>
<UL>
<DL>
<DT><b>--in</b> &#60;<i>input source and target</i>&#62;
<DD> These two strings specify the names of the source and target textures for optical flow.<br>
The images should have the same resolution and should be represented in either PNG or JPEG format.<br>
Alternatively, if colors are stored directly per-vertex, the two strings should give the names of the files storing the surface geometry, in <A HREF="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format, with the files representing the same geometry but with different colors at the vertices.

<DT>[<b>--mesh</b> &#60;<i>base geometry</i>&#62;]
<DD> This string specifies the name of the file storing the surface geometry and associated textured coordinates.<BR>
The geometry is assumed to be in <A HREF="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format, with texture coordinates stored either at the vertices or at the triangle corners.<br>

<DT>[<b>--levels</b> &#60;<i>smoothing levels</i>&#62;]
<DD> This integer gives the number of (multi-resolution) hierarchy levels over which the optical flow should be run.<br>
The default value for this parameter is 7.<br>

<DT>[<b>--sSmooth</b> &#60;<i>smoothing weight</i>&#62;]
<DD> This floating point value gives the smoothing weight for scalar fields at the coarsest resolution of the hierarchy.<br>
The default value for this parameter is 0.003.<br>

<DT>[<b>--vfSmooth</b> &#60;<i>smoothing weight</i>&#62;]
<DD> This floating point value gives the smoothing weight for vector fields at the coarsest resolution of the hierarchy.<br>
The default value for this parameter is 1000000.<br>

<DT>[<b>--out</b> &#60;<i>output image headers</i>&#62;]
<DD> This string gives the header for the files into which the in-between textures are written.<br>

<DT>[<b>--ext</b> &#60;<i>output image extension</i>&#62;]
<DD> This string gives the extension for the files into which the in-between textures are written. (Only PNG and JPEG are supported.)<br>
The default value for this parameter is "png".<br> 

<DT>[<b>--keyFrames</b> &#60;<i>number of textures to output</i>&#62;]
<DD> This integer gives the number of in-between textures that should be generated.<br>
The default value for this parameter is 2.<br>

<DT>[<b>--padRadius</b> &#60;<i>texture padding radius</i>&#62;]
<DD> This integer gives the dilation radius (in texels) around the rasterization of the mesh into texture space for which texture values should be computed.<br>
The default value for this parameter is 3.<br>

<DT>[<b>--eLength</b> &#60;<i>target edge length</i>&#62;]
<DD> This floating point value specifies the target edge length (as a fraction of the total circumference) up to which the surface should be subdivided for sampling the texture to the vertices.<br>
The default value for this parameter is 0.008.<br>

<DT>[<b>--dogWeight</b> &#60;<i>difference of Gaussians blending weight</i>&#62;]
<DD> This floating point value gives the relative importance of matching the difference of Gaussians of pixel values to matching the actual pixel values when estimating the optical flow field.<br>
The default value for this parameter is 1.0.<br>
  
<DT>[<b>--whitney</b>]
<DD> Enabling this flag will use the Whitney (1-form) basis functions for representing the flow field instead of the conforming basis.<br>

<DT>[<b>--error</b>]
<DD> Enabling this flag output the error in the fit between advected textures.<br>

<DT>[<b>--verbose</b>]
<DD> Enabling this flag provides a verbose description of the break-down in running times.<br>

</UL>

<HR>
<A NAME="CHANGES"><B>CHANGES</B></A><br>
<!--
<A HREF="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version3/">Version 3</A>:
<OL>
<LI> The implementation of the <b>--samplesPerNode</b> parameter has been modified so that a value of "1" more closely corresponds to a distribution with one sample per leaf node.
<LI> The code has been modified to support compilation under MSVC 2010 and the associated solution and project files are now provided. (Due to a bug in the Visual Studios compiler, this required modifying the implementation of some of the bit-shifting operators.)
</OL>
-->
