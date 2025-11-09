#ifndef __Elasticity_Custom_Boundary_h__
#define __Elasticity_Custom_Boundary_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "ElasticityBase.h"
#include "Elasticity_Kugelstadt2021.h"
#include "SPlisHSPlasH/Utilities/MatrixFreeSolver.h"
#if USE_AVX
#include "SPlisHSPlasH/Utilities/AVX_math.h"
#include "SPlisHSPlasH/Utilities/CholeskyAVXSolver.h"
#endif


namespace SPH
{
	/** \brief This class inherits the implicit SPH formulation for
	* incompressible linearly elastic solids introduced
	* by Kugelstadt et al. [KBF+21] and loads custom boundary conditions.
	*
	* References:
	* - [KBF+21] Tassilo Kugelstadt, Jan Bender, José Antonio Fernández-Fernández, Stefan Rhys Jeske, Fabian Löschner, Andreas Longva. 
	* Fast Corotated Elastic SPH Solids with Implicit Zero-Energy Mode Control. Proceedings of the ACM on Computer Graphics and 
	* Interactive Techniques, 2021. URL: http://dx.doi.org/10.1145/3480142
	*/
	class Elasticity_Custom_Boundary : public Elasticity_Kugelstadt2021
	{
	protected:

		virtual void initParameters();
		virtual void determineFixedParticles();
		void updateBoundaryParticles();
		
		// Cached boundary animation data
		bool m_bdy_loaded = false;
		unsigned int m_bdy_frame_index = 0u;
		unsigned int m_bdy_num_frames = 0u;
		unsigned int m_bdy_num_boundary = 0u;
		std::vector<unsigned int> m_bdy_indices;
		std::vector<Vector3r> m_bdy_pos_seq;		// flattened as [frame * num_boundary + k]

	public:

        static int BOUNDARY_DATA_FILE;

        std::string m_boundary_data_file; 

		Elasticity_Custom_Boundary(FluidModel *model);
		virtual ~Elasticity_Custom_Boundary(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new Elasticity_Custom_Boundary(model); }

		virtual void step();
	};
}

#endif
