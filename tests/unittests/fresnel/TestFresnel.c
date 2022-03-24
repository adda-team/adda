#include <check.h>
#include <stdio.h>
#include "../../../src/cmplx.h"

#define ck_assert_doublecomplex_eq_tol(a, b, tol)                                         \
	do {                                                                                  \
		ck_assert_double_eq_tol(creal((doublecomplex) a), creal((doublecomplex) b), tol); \
		ck_assert_double_eq_tol(cimag((doublecomplex) a), cimag((doublecomplex) b), tol); \
	} while (false)

void PrintError(const char * restrict fmt, ... ) {
	printf(fmt);
};

START_TEST(test_FresnelS) {
	double tol = 1e-16;
	doublecomplex t = FresnelTS(1.3 + I * 0.1, 1.3 + I * 0.1);
	ck_assert_doublecomplex_eq_tol(t, 1, tol);
	doublecomplex r = FresnelRS(1.3 + I * 0.1, 1.3 + I * 0.1);
	ck_assert_doublecomplex_eq_tol(r, 0, tol);

}
END_TEST

START_TEST(test_FresnelP) {
	double tol = 1e-16;
	doublecomplex t = FresnelTP(1.3 + I * 0.1, 1.3 + I * 0.1, 1);
	ck_assert_doublecomplex_eq_tol(t, 1, tol);
	doublecomplex r = FresnelRP(1.3 + I * 0.1, 1.3 + I * 0.1, 1);
	ck_assert_doublecomplex_eq_tol(t, 1, tol);

}
END_TEST

void check_all_cases_SubstrateFresnel (const struct Substrate sub, const double wave_num, const bool z_positive,
		                               const doublecomplex sqr_long_k, const doublecomplex ki, 
									   const doublecomplex ts_ref, const doublecomplex rs_ref, 
		                               const doublecomplex tp_ref, const doublecomplex rp_ref, doublecomplex kt_ref,
		                               const double tol) {
#define ARG_SEQ_FRESNEL sub, wave_num, z_positive, sqr_long_k, ki
#define CHECK_FRESNEL_HELPER(c) ck_assert_doublecomplex_eq_tol(c, c##_ref, tol)
#define REINIT_COEF ts = 2, rs = 2, tp = 2, rp = 2, kt = 2
	
	doublecomplex ts, rs, tp, rp, kt;
	REINIT_COEF;
	
	SubstrateFresnel(ARG_SEQ_FRESNEL, &ts, NULL, NULL, NULL, &kt);
	CHECK_FRESNEL_HELPER(ts);
	CHECK_FRESNEL_HELPER(kt);
	REINIT_COEF;
	
	SubstrateFresnel(ARG_SEQ_FRESNEL, NULL, &rs, NULL, NULL, NULL);
	CHECK_FRESNEL_HELPER(rs);
	REINIT_COEF;
	SubstrateFresnel(ARG_SEQ_FRESNEL, NULL, NULL, &tp, NULL, NULL);
	CHECK_FRESNEL_HELPER(tp);
	REINIT_COEF;

	SubstrateFresnel(ARG_SEQ_FRESNEL, NULL, NULL, NULL, &rp, NULL);
	CHECK_FRESNEL_HELPER(rp);
	REINIT_COEF;

	SubstrateFresnel(ARG_SEQ_FRESNEL, &ts, &rs, NULL, NULL, NULL);
	CHECK_FRESNEL_HELPER(ts);
	CHECK_FRESNEL_HELPER(rs);
	REINIT_COEF;

	SubstrateFresnel(ARG_SEQ_FRESNEL, NULL, NULL, &tp, &rp, NULL);
	CHECK_FRESNEL_HELPER(tp);
	CHECK_FRESNEL_HELPER(rp);
	REINIT_COEF;

	SubstrateFresnel(ARG_SEQ_FRESNEL, &ts, &rs, &tp, &rp, &kt);
	CHECK_FRESNEL_HELPER(ts);
	CHECK_FRESNEL_HELPER(rs);
	CHECK_FRESNEL_HELPER(tp);
	CHECK_FRESNEL_HELPER(rp);
	CHECK_FRESNEL_HELPER(kt);
	REINIT_COEF;
	
	
#undef ARG_SEQ_FRESNEL
#undef CHECK_FRESNEL_HELPER
#undef REINIT_COEF
}

START_TEST(test_SubstrateFresnel_1) {
	struct Substrate sub;
	sub.mInf = false;
	sub.m[0] = 1.3;
	sub.m[1] = 1.3;
	sub.h[0] = 2;
	sub.N = 2;
	doublecomplex ki = 1;
	doublecomplex kt_ref = CalculateKt(ki, 1, sub.m[1], 0);
	doublecomplex ts_ref = FresnelTS(ki, kt_ref);
	doublecomplex tp_ref = FresnelTP(ki, kt_ref, 1.3);
	doublecomplex rs_ref = FresnelRS(ki, kt_ref);
	doublecomplex rp_ref = FresnelRP(ki, kt_ref, 1.3);
	check_all_cases_SubstrateFresnel(sub, 1, false, 0, ki,
                                     ts_ref, rs_ref, tp_ref, rp_ref, kt_ref, 1e-14);

}
END_TEST

START_TEST(test_SubstrateFresnel_2) {
    struct Substrate sub;
    sub.mInf = false;
    sub.m[0] = 1;
    sub.m[1] = 1.3;
    sub.h[0] = 2;
    sub.N = 2;
    doublecomplex ki = 1;
    double wave_num = 1;
    doublecomplex eL_sqr = cexp(I * wave_num * sub.h[0] * ki * 2);
    doublecomplex kt_ref = CalculateKt(ki, 1, sub.m[1], 0);
    doublecomplex ts_ref = FresnelTS(ki, kt_ref);
    doublecomplex tp_ref = FresnelTP(ki, kt_ref, 1.3);
    doublecomplex rs_ref = FresnelRS(ki, kt_ref) * eL_sqr;
    doublecomplex rp_ref = FresnelRP(ki, kt_ref, 1.3) * eL_sqr;
    check_all_cases_SubstrateFresnel(sub, wave_num, false, 0, ki,
                                     ts_ref, rs_ref, tp_ref, rp_ref, kt_ref, 1e-14);

}
END_TEST

START_TEST(test_SubstrateFresnel_3) {
    struct Substrate sub;
    sub.mInf = false;
    sub.m[0] = 1.3;
    sub.m[1] = 1.3;
    sub.h[0] = 2;
    sub.N = 2;
    doublecomplex ki = 1.3;
    double wave_num = 1;
    doublecomplex eL_sqr = cexp(I * wave_num * sub.h[0] * ki * 2);
    doublecomplex kt_ref = CalculateKt(ki, 1.3, 1, 0);
    doublecomplex ts_ref = FresnelTS(ki, kt_ref);
    doublecomplex tp_ref = FresnelTP(ki, kt_ref, 1/1.3);
    doublecomplex rs_ref = FresnelRS(ki, kt_ref) * eL_sqr;
    doublecomplex rp_ref = FresnelRP(ki, kt_ref, 1/1.3) * eL_sqr;
    check_all_cases_SubstrateFresnel(sub, 1, true, 0, ki,
                                     ts_ref, rs_ref, tp_ref, rp_ref, kt_ref, 1e-14);

}
END_TEST

Suite * fresnel_suite(void)
{
	Suite *s;
	TCase *tc_Fresnel;

	s = suite_create("Fresnel");

	tc_Fresnel = tcase_create("Fresnel");


	tcase_add_test(tc_Fresnel, test_FresnelS);
	tcase_add_test(tc_Fresnel, test_FresnelP);
	tcase_add_test(tc_Fresnel, test_SubstrateFresnel_1);
    tcase_add_test(tc_Fresnel, test_SubstrateFresnel_2);
    tcase_add_test(tc_Fresnel, test_SubstrateFresnel_3);
	suite_add_tcase(s, tc_Fresnel);

	return s;
}

int main(void)
{
	int number_failed;
	Suite *s;
	SRunner *sr;

	s = fresnel_suite();
	sr = srunner_create(s);

	srunner_run_all(sr, CK_NORMAL);
	number_failed = srunner_ntests_failed(sr);
	srunner_free(sr);
	return (number_failed == 0) ? 0 : 1;
}