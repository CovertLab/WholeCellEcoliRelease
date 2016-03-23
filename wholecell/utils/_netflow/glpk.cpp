#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python.hpp>
#include <boost/numpy.hpp>
#include <iostream>
#include <glpk.h>

using namespace boost::python;
using namespace boost::numpy;

#define ASSERT_THROW(a, msg) if (!(a)) throw std::runtime_error(msg)

class glpk
{
	public:
		glpk()
		{
			this->lp = glp_create_prob();
			glp_init_smcp(&(this->lp_params));
			this->n_eq_constraints = 0;
			this->n_vars = 0;
		}
		~glpk(){
			glp_delete_prob(this->lp);
		}
		void add_rows(int n_rows)
		{
			glp_add_rows(this->lp, n_rows);
			this->n_eq_constraints += n_rows;
		}
		void add_cols(int n_cols)
		{
			glp_add_cols(this->lp, n_cols);
			this->n_vars += n_cols;
		}
		void set_col_bounds(int index, double lower, double upper)
		{
			if(lower == upper){
				glp_set_col_bnds(this->lp, index, GLP_FX, lower, upper);
				return;
			}
			if(!isinf(lower) && !isinf(upper)){
				glp_set_col_bnds(this->lp, index, GLP_DB, lower, upper);
				return;
			}
			if(isinf(lower) && isinf(upper)){
				glp_set_col_bnds(this->lp, index, GLP_FR, lower, upper);
				return;
			}
			if(isinf(upper)){
				glp_set_col_bnds(this->lp, index, GLP_LO, lower, upper);
				return;
			}
			if(isinf(lower)){
				glp_set_col_bnds(this->lp, index, GLP_UP, lower, upper);
				return;
			}
		}
		void set_obj_coef(int index, double value)
		{
			glp_set_obj_coef(this->lp, index, value);
		}
		void set_sense_max()
		{
			glp_set_obj_dir(lp, GLP_MAX);
		}
		void set_sense_min()
		{
			glp_set_obj_dir(lp, GLP_MIN);
		}
		void set_quiet()
		{
			this->lp_params.msg_lev = GLP_MSG_ERR;
		}
		void set_quiet_quiet()
		{
			this->lp_params.msg_lev = GLP_MSG_OFF;
		}
		void set_verbose()
		{
			this->lp_params.msg_lev = GLP_MSG_ON;
		}
		void set_verbose_verbose()
		{
			this->lp_params.msg_lev = GLP_MSG_ALL;
		}
		void optimize()
		{
        		int status = glp_simplex(this->lp, &(this->lp_params));
			switch(status){
				case 0: break;
				case GLP_EBADB: throw std::runtime_error("GLP_EBADB"); break;
				case GLP_ESING: throw std::runtime_error("GLP_ESING"); break;
				case GLP_ECOND: throw std::runtime_error("GLP_ECOND"); break;
				case GLP_EBOUND: throw std::runtime_error("GLP_EBOUND"); break;
				case GLP_EFAIL: throw std::runtime_error("GLP_EFAIL"); break;
				case GLP_EOBJLL: throw std::runtime_error("GLP_EOBJLL"); break;
				case GLP_EOBJUL: throw std::runtime_error("GLP_EOBJUL"); break;
				case GLP_EITLIM: throw std::runtime_error("GLP_EITLIM"); break;
				case GLP_ETMLIM: throw std::runtime_error("GLP_ETMLIM"); break;
				case GLP_ENOPFS: throw std::runtime_error("GLP_ENOPFS"); break;
				case GLP_ENODFS: throw std::runtime_error("GLP_ENODFS"); break;
				default: throw std::runtime_error("UNKNOWN SOLVER RETURN VALUE"); break;
			}
			if(glp_get_status(this->lp) != GLP_OPT)
				throw std::runtime_error("NON-OPTIMAL SOLUTION");
		}
		double get_primal_value(int index)
		{
			return glp_get_col_prim(this->lp, index);
		}
		double get_objective_value()
		{
			return glp_get_obj_val(this->lp);
		}
		void solution(ndarray& objective_value_array, ndarray& primal_variable_array)
		{
			ASSERT_THROW((objective_value_array.get_dtype() == dtype::get_builtin<double>()), "Expected array of type np.float64");
			ASSERT_THROW((primal_variable_array.get_dtype() == dtype::get_builtin<double>()), "Expected array of type np.float64");
			ASSERT_THROW((objective_value_array.shape(0) == 1), "Must be of size 1!");
			ASSERT_THROW((primal_variable_array.shape(0) == this->n_vars + 1), "Must be of size n_vars + 1 (due to zero vs one indexing)!");
			double* objective_value = (double *) objective_value_array.get_data();
			double* primal_variable = (double *) primal_variable_array.get_data();
        		int status = glp_simplex(this->lp, &(this->lp_params));
			switch(status){
				case 0: break;
				case GLP_EBADB: throw std::runtime_error("GLP_EBADB"); break;
				case GLP_ESING: throw std::runtime_error("GLP_ESING"); break;
				case GLP_ECOND: throw std::runtime_error("GLP_ECOND"); break;
				case GLP_EBOUND: throw std::runtime_error("GLP_EBOUND"); break;
				case GLP_EFAIL: throw std::runtime_error("GLP_EFAIL"); break;
				case GLP_EOBJLL: throw std::runtime_error("GLP_EOBJLL"); break;
				case GLP_EOBJUL: throw std::runtime_error("GLP_EOBJUL"); break;
				case GLP_EITLIM: throw std::runtime_error("GLP_EITLIM"); break;
				case GLP_ETMLIM: throw std::runtime_error("GLP_ETMLIM"); break;
				case GLP_ENOPFS: throw std::runtime_error("GLP_ENOPFS"); break;
				case GLP_ENODFS: throw std::runtime_error("GLP_ENODFS"); break;
				default: throw std::runtime_error("UNKNOWN SOLVER RETURN VALUE"); break;
			}
			if(glp_get_status(this->lp) != GLP_OPT)
				throw std::runtime_error("NON-OPTIMAL SOLUTION");
			objective_value[0] = glp_get_obj_val(this->lp);
			for(int i = 1; i <= n_vars; i++)
				primal_variable[i] = glp_get_col_prim(this->lp, i);
		}
		void add_eq_constrs(ndarray& ia_array, ndarray& ja_array, ndarray& ar_array)
		{
			ASSERT_THROW((ia_array.get_dtype() == dtype::get_builtin<int32_t>()), "Expected array of type np.int32");
			ASSERT_THROW((ja_array.get_dtype() == dtype::get_builtin<int32_t>()), "Expected array of type np.int32");
			ASSERT_THROW((ar_array.get_dtype() == dtype::get_builtin<double>()), "Expected array of type np.float64");
			ASSERT_THROW((ia_array.shape(0) == ja_array.shape(0)), "Sizes must match!");
			ASSERT_THROW((ia_array.shape(0) == ar_array.shape(0)), "Sizes must match!");
			int32_t* ia = (int32_t *) ia_array.get_data();
			int32_t* ja = (int32_t *) ja_array.get_data();
			double* ar = (double *) ar_array.get_data();
			for(int row = 1; row <= this->n_eq_constraints; row++)
				glp_set_row_bnds(this->lp, row, GLP_FX, 0.0, 0.0);
			int32_t n_elems = ia_array.shape(0);
			glp_load_matrix(this->lp, n_elems - 1, ia, ja, ar);
			

		}
		void set_mat_row(int row, int len, ndarray& ja_array, ndarray& ar_array)
		{
			ASSERT_THROW((ja_array.get_dtype() == dtype::get_builtin<int32_t>()), "Expected array of type np.int32");
			ASSERT_THROW((ar_array.get_dtype() == dtype::get_builtin<double>()), "Expected array of type np.float64");
			ASSERT_THROW((ja_array.shape(0) == ar_array.shape(0)), "Sizes must match!");
			int32_t* ja = (int32_t *) ja_array.get_data();
			double* ar = (double *) ar_array.get_data();

			glp_set_mat_row(this->lp, row, len, ja, ar);
		}
	private:
		glp_prob *lp;
		glp_smcp lp_params;
		int n_eq_constraints;
		int n_vars;
};

BOOST_PYTHON_MODULE(glpk)
{
    Py_Initialize();
    initialize();

    class_<glpk>("glpk")
	.def("add_rows", &glpk::add_rows)
	.def("add_cols", &glpk::add_cols)
	.def("set_col_bounds", &glpk::set_col_bounds)
	.def("set_obj_coef", &glpk::set_obj_coef)
	.def("add_eq_constrs", &glpk::add_eq_constrs)
	.def("set_mat_row", &glpk::set_mat_row)
	.def("set_sense_max", &glpk::set_sense_max)
	.def("set_sense_min", &glpk::set_sense_min)
	.def("set_quiet", &glpk::set_quiet)
	.def("set_quiet_quiet", &glpk::set_quiet_quiet)
	.def("set_verbose", &glpk::set_verbose)
	.def("set_verbose_verbose", &glpk::set_verbose_verbose)
	.def("optimize", &glpk::optimize)
	.def("get_primal_value", &glpk::get_primal_value)
	.def("get_objective_value", &glpk::get_objective_value)
	.def("solution", &glpk::solution);
}

