#define PY_SSIZE_T_CLEAN
#include <Python.h>
#ifndef OVERRIDE_STDC_TEST
#define OVERRIDE_STDC_TEST
#endif
#include "../../../src/cmplx.h"

#define FRESNEL_DECL(POLAR_TYPE) {            \
	"Fresnel"#POLAR_TYPE,                     \
	Py_Fresnel##POLAR_TYPE,                   \
	METH_VARARGS,                             \
	"Python interface for Fresnel"#POLAR_TYPE \
	}

#define FRESNEL_ARGNUM_S 2
#define FRESNEL_ARGNUM_P 3
#define FRESNEL_ARRAY_SEQ_S(a) a[0], a[1]
#define FRESNEL_ARRAY_SEQ_P(a) a[0], a[1], a[2]
#define FRESNEL_FORMAT_S "DD"
#define FRESNEL_FORMAT_P "DDD"

#define MAKE_FRESNEL(T_OR_R, S_OR_P)                                                     \
    static PyObject * Py_Fresnel##T_OR_R##S_OR_P(PyObject * self, PyObject * args) {     \
		doublecomplex c_arg[FRESNEL_ARGNUM_##S_OR_P];                                    \
		Py_complex py_arg[FRESNEL_ARGNUM_##S_OR_P];                                      \
		if (!PyArg_ParseTuple(args, FRESNEL_FORMAT_##S_OR_P,                             \
			                  FRESNEL_ARRAY_SEQ_##S_OR_P(&py_arg))                       \
		)                                                                                \
			return NULL;                                                                 \
		for (int i = 0; i < FRESNEL_ARGNUM_##S_OR_P; ++i)                                \
			c_arg[i] = py_arg[i].real + I * py_arg[i].imag;                              \
        doublecomplex res = Fresnel##T_OR_R##S_OR_P(FRESNEL_ARRAY_SEQ_##S_OR_P(c_arg));  \
		Py_complex py_res = {creal(res), cimag(res)};                                    \
		return PyComplex_FromCComplex(py_res);                                           \
	}

MAKE_FRESNEL(T, S);
MAKE_FRESNEL(R, S);
MAKE_FRESNEL(T, P);
MAKE_FRESNEL(R, P);

void PrintError(const char * restrict fmt, ... ) {
	printf(fmt);
};

#TODO THIS WRAPPER IS INCOMPLETE
static PyObject * Py_SubstrateFresnel(PyObject * self, PyObject * args) {
	PyObject  py_sub;
	struct Substrate c_sub;
	Py_complex py_wave_num;
	doublecomplex c_wave_num;
	bool is_z_positive;
	Py_complex py_sqr_long_k;
	doublecomplex c_sqr_long_k;
	Py_complex py_ki;
	doublecomplex c_ki;
	PyObject py_coef_set;
	doublecomplex c_coef[4];
	doublecomplex * c_coef_ptr[4];
	for (int i = 0; i < 4; ++i)
		c_coef_ptr[i] = NULL;

	if (!PyArg_ParseTuple(args, "ODpDDO", &py_sub, &py_wave_num, &is_z_positive, &py_sqr_long_k, &py_ki, &py_coef_set))
		return NULL;
	c_wave_num = py_wave_num.real + I * py_wave_num.imag;
	c_sqr_long_k = py_sqr_long_k.real + I * py_sqr_long_k.imag;
	c_ki = py_ki.real + I * py_ki.imag;
	PyObject * tmp_obj;
	tmp_obj = PyDict_GetItemString(&py_sub, "mInf");
	if (!PyArg_ParseTuple(tmp_obj, "p", &c_sub.mInf))
		return NULL;
	tmp_obj = PyDict_GetItemString(&py_sub, "N");
	if (!PyArg_ParseTuple(tmp_obj, "i", &c_sub.N))
		return NULL;
	tmp_obj = PyDict_GetItemString(&py_sub, "hP");
	if (!PyArg_ParseTuple(tmp_obj, "d", &c_sub.hP))
		return NULL;
	Py_complex m[2];
	tmp_obj = PyDict_GetItemString(&py_sub, "m");
	if (!PyArg_ParseTuple(tmp_obj, "(DD)", m))
		return NULL;
	c_sub.m[0] = m[0].real + I * m[0].imag;
	c_sub.m[1] = m[0].real + I * m[1].imag;
	tmp_obj = PyDict_GetItemString(&py_sub, "h");
	if (!PyArg_ParseTuple(tmp_obj, "(dd)", c_sub.h))
		return NULL;
	SubstrateFresnel(c_sub, c_wave_num, is_z_positive, c_sqr_long_k, c_ki,
					 &c_coef[0], &c_coef[1], &c_coef[2], &c_coef[3], NULL);
	Py_complex py_res = {creal(c_coef[0]), cimag(c_coef[0])};
		return PyComplex_FromCComplex(py_res);
}

static PyMethodDef ADDAWrappersMethods[] = {
		FRESNEL_DECL(TS),
		FRESNEL_DECL(RS),
		FRESNEL_DECL(TP),
		FRESNEL_DECL(RP),
		{"SubstrateFresnel", Py_SubstrateFresnel, METH_VARARGS, "Python interface for SubstrateFresnel"},
		{NULL, NULL, 0, NULL}
};

static struct PyModuleDef ADDAWrappersModule = {
		PyModuleDef_HEAD_INIT,
		"Fresnel",
		"Python interface for fresnel coefficients",
		-1,
		ADDAWrappersMethods
};

PyMODINIT_FUNC PyInit_ADDAWrappers(void) {
	return PyModule_Create(&ADDAWrappersModule);
}

