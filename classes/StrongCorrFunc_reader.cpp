#include "StrongCorrFunc_reader.h"

StrongCorrFunc_reader::StrongCorrFunc_reader(const char* filename)
{
  CC1_array = NULL;
  CC2_array = NULL;
  CC3_array = NULL;
  CC4_real_array = NULL;
  CC4_imag_array = NULL;
  InitTable(filename);
}

StrongCorrFunc_reader::~StrongCorrFunc_reader()
{
  if(CC1_array) delete [] CC1_array; 
  if(CC2_array) delete [] CC2_array;
  if(CC3_array) delete [] CC3_array;
  if(CC4_real_array) delete [] CC4_real_array;
  if(CC4_imag_array) delete [] CC4_imag_array;
}

double StrongCorrFunc_reader::getValue(const double alpha, const double k, const double Rcc) const
{
  return get_value_CC1(alpha, k, Rcc) + get_value_CC2(alpha, k, Rcc);
}

double StrongCorrFunc_reader::read_table1(const double alpha, const double k, const double Rcc) const
{
  if(!CC1_array)
    return 0;
  int ialpha = floor((alpha - alpha_min) / d_alpha); 
  int ik = floor((k - k_min) / d_k); 
  int iRcc = floor((Rcc - Rcc_min) / d_Rcc); 
  if(ialpha >= 0 && ialpha < (Nalpha + 1))
    if(ik >= 0 && ik < (Nk + 1))
      if(iRcc >= 0 && iRcc < (NRcc + 1))
      return CC1_array[ialpha][ik][iRcc];
  return 0;
}

double StrongCorrFunc_reader::read_table2(const double alpha, const double k, const double Rcc) const
{
  if(!CC2_array)
    return 0;
  int ialpha = floor((alpha - alpha_min) / d_alpha); 
  int ik = floor((k - k_min) / d_k); 
  int iRcc = floor((Rcc - Rcc_min) / d_Rcc); 
  if(ialpha >= 0 && ialpha < (Nalpha + 1))
    if(ik >= 0 && ik < (Nk + 1))
      if(iRcc >= 0 && iRcc < (NRcc + 1))
      return CC2_array[ialpha][ik][iRcc];
  return 0;
}

double StrongCorrFunc_reader::get_value_CC1(const double alpha, const double k, const double Rcc) const
{
  if(!CC1_array)
    return 0;
  int ik = floor((k - k_min) / d_k);
  double klo = k_min + ik * d_k;
  int ialpha = floor((alpha - alpha_min) / d_alpha);
  double alphalo = alpha_min + ialpha * d_alpha;
  int iRcc = floor((Rcc - Rcc_min) / d_Rcc);
  double Rcclo = Rcc_min + iRcc * d_Rcc;
  if(ik < 0 || ik > Nk || ialpha < 0 || ialpha > Nalpha || iRcc < 0 || iRcc > NRcc)
    return 0;
  if(ik < Nk && ialpha < Nalpha && iRcc < NRcc)
  { 
    double klo_lint_iRcc = CC1_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC1_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
    double khi_lint_iRcc = CC1_array[ialpha + 1][ik][iRcc] * (klo - k + d_k) / d_k + CC1_array[ialpha + 1][ik + 1][iRcc] * (k - klo) /d_k;
    double klo_lint_iRcc_plus = CC1_array[ialpha][ik][iRcc + 1] * (klo - k + d_k) / d_k + CC1_array[ialpha][ik + 1][iRcc + 1] * (k - klo) / d_k;
    double khi_lint_iRcc_plus = CC1_array[ialpha + 1][ik][iRcc + 1] * (klo - k + d_k) / d_k + CC1_array[ialpha + 1][ik + 1][iRcc + 1] * (k - klo) /d_k;

    double alphalo_lint = klo_lint_iRcc * (alphalo - alpha + d_alpha) / d_alpha + khi_lint_iRcc * (alpha - alphalo) / d_alpha;
    double alphahi_lint = klo_lint_iRcc_plus * (alphalo - alpha + d_alpha) / d_alpha + khi_lint_iRcc_plus * (alpha - alphalo) / d_alpha;
    return alphalo_lint * (Rcclo - Rcc + d_Rcc) / d_Rcc + alphahi_lint * (Rcc - Rcclo) / d_Rcc;
  }
  if(ik == Nk && ialpha < Nalpha && iRcc < NRcc)
  {
    double Rcclo_lint = CC1_array[ialpha][ik][iRcc] * (Rcclo - Rcc + d_Rcc) / d_Rcc + CC1_array[ialpha][ik][iRcc + 1] * (Rcc - Rcclo) / d_Rcc;
    double Rcchi_lint = CC1_array[ialpha + 1][ik][iRcc] * (Rcclo - Rcc + d_Rcc) / d_Rcc + CC1_array[ialpha + 1][ik][iRcc + 1] * (Rcc - Rcclo) /d_Rcc;
    return Rcclo_lint * (alphalo - alpha + d_alpha) / d_alpha + Rcchi_lint * (alpha - alphalo) / d_alpha;
  }
  if(ialpha == Nalpha && ik < Nk && iRcc < NRcc)
  {
    double klo_lint = CC1_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC1_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
    double khi_lint = CC1_array[ialpha][ik][iRcc + 1] * (klo - k + d_k) / d_k + CC1_array[ialpha][ik + 1][iRcc + 1] * (k - klo) /d_k;
    return klo_lint * (Rcclo - Rcc + d_Rcc) / d_Rcc + khi_lint * (Rcc - Rcclo) / d_Rcc;
  }
  if(iRcc == NRcc && ik < Nk && ialpha < Nalpha)
  {
    double klo_lint = CC1_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC1_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
    double khi_lint = CC1_array[ialpha + 1][ik][iRcc] * (klo - k + d_k) / d_k + CC1_array[ialpha + 1][ik + 1][iRcc] * (k - klo) /d_k;
    return klo_lint * (alphalo - alpha + d_alpha) / d_alpha + khi_lint * (alpha - alphalo) / d_alpha;
  }
  else if(iRcc == NRcc && ik == Nk && ialpha < Nalpha)
    return CC1_array[ialpha][ik][iRcc] * (alphalo - alpha + d_alpha) / d_alpha + CC1_array[ialpha + 1][ik][iRcc] * (alpha - alphalo) / d_alpha;
  else if(iRcc == NRcc && ik < Nk && ialpha == Nalpha)
    return CC1_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC1_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
  else if(iRcc < NRcc && ik == Nk && ialpha == Nalpha)
    return CC1_array[ialpha][ik][iRcc] * (Rcclo - Rcc + d_Rcc) / d_Rcc + CC1_array[ialpha][ik][iRcc + 1] * (Rcc - Rcclo) / d_Rcc;
  else if(ik == Nk && ialpha == Nalpha && iRcc == NRcc)
    return CC1_array[ialpha][ik][iRcc];
  return 0;
}

double StrongCorrFunc_reader::get_value_CC2(const double alpha, const double k, const double Rcc) const
{
  if(!CC2_array)
    return 0;
  int ik = floor((k - k_min) / d_k);
  double klo = k_min + ik * d_k;
  int ialpha = floor((alpha - alpha_min) / d_alpha);
  double alphalo = alpha_min + ialpha * d_alpha;
  int iRcc = floor((Rcc - Rcc_min) / d_Rcc);
  double Rcclo = Rcc_min + iRcc * d_Rcc;
  if(ik < 0 || ik > Nk || ialpha < 0 || ialpha > Nalpha || iRcc < 0 || iRcc > NRcc)
    return 0;
  if(ik < Nk && ialpha < Nalpha && iRcc < NRcc)
  { 
    double klo_lint_iRcc = CC2_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC2_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
    double khi_lint_iRcc = CC2_array[ialpha + 1][ik][iRcc] * (klo - k + d_k) / d_k + CC2_array[ialpha + 1][ik + 1][iRcc] * (k - klo) /d_k;
    double klo_lint_iRcc_plus = CC2_array[ialpha][ik][iRcc + 1] * (klo - k + d_k) / d_k + CC2_array[ialpha][ik + 1][iRcc + 1] * (k - klo) / d_k;
    double khi_lint_iRcc_plus = CC2_array[ialpha + 1][ik][iRcc + 1] * (klo - k + d_k) / d_k + CC2_array[ialpha + 1][ik + 1][iRcc + 1] * (k - klo) /d_k;

    double alphalo_lint = klo_lint_iRcc * (alphalo - alpha + d_alpha) / d_alpha + khi_lint_iRcc * (alpha - alphalo) / d_alpha;
    double alphahi_lint = klo_lint_iRcc_plus * (alphalo - alpha + d_alpha) / d_alpha + khi_lint_iRcc_plus * (alpha - alphalo) / d_alpha;
    return alphalo_lint * (Rcclo - Rcc + d_Rcc) / d_Rcc + alphahi_lint * (Rcc - Rcclo) / d_Rcc;
  }
  else if(ik == Nk && ialpha < Nalpha && iRcc < NRcc)
  {
    double Rcclo_lint = CC2_array[ialpha][ik][iRcc] * (Rcclo - Rcc + d_Rcc) / d_Rcc + CC2_array[ialpha][ik][iRcc + 1] * (Rcc - Rcclo) / d_Rcc;
    double Rcchi_lint = CC2_array[ialpha + 1][ik][iRcc] * (Rcclo - Rcc + d_Rcc) / d_Rcc + CC2_array[ialpha + 1][ik][iRcc + 1] * (Rcc - Rcclo) /d_Rcc;
    return Rcclo_lint * (alphalo - alpha + d_alpha) / d_alpha + Rcchi_lint * (alpha - alphalo) / d_alpha;
  }
  else if(ialpha == Nalpha && ik < Nk && iRcc < NRcc)
  {
    double klo_lint = CC2_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC2_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
    double khi_lint = CC2_array[ialpha][ik][iRcc + 1] * (klo - k + d_k) / d_k + CC2_array[ialpha][ik + 1][iRcc + 1] * (k - klo) /d_k;
    return klo_lint * (Rcclo - Rcc + d_Rcc) / d_Rcc + khi_lint * (Rcc - Rcclo) / d_Rcc;
  }
  else if(iRcc == NRcc && ik < Nk && ialpha < Nalpha)
  {
    double klo_lint = CC2_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC2_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
    double khi_lint = CC2_array[ialpha + 1][ik][iRcc] * (klo - k + d_k) / d_k + CC2_array[ialpha + 1][ik + 1][iRcc] * (k - klo) /d_k;
    return klo_lint * (alphalo - alpha + d_alpha) / d_alpha + khi_lint * (alpha - alphalo) / d_alpha;
  }
  else if(iRcc == NRcc && ik == Nk && ialpha < Nalpha)
    return CC2_array[ialpha][ik][iRcc] * (alphalo - alpha + d_alpha) / d_alpha + CC2_array[ialpha + 1][ik][iRcc] * (alpha - alphalo) / d_alpha;
  else if(iRcc == NRcc && ik < Nk && ialpha == Nalpha)
    return CC2_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC2_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
  else if(iRcc < NRcc && ik == Nk && ialpha == Nalpha)
    return CC2_array[ialpha][ik][iRcc] * (Rcclo - Rcc + d_Rcc) / d_Rcc + CC2_array[ialpha][ik][iRcc + 1] * (Rcc - Rcclo) / d_Rcc;
  else if(ik == Nk && ialpha == Nalpha && iRcc == NRcc)
    return CC2_array[ialpha][ik][iRcc];
  return 0;
}

double StrongCorrFunc_reader::get_value_CC3(const double alpha, const double k, const double Rcc) const
{
  if(!CC3_array)
    return 0;
  int ik = floor((k - k_min) / d_k);
  double klo = k_min + ik * d_k;
  int ialpha = floor((alpha - alpha_min) / d_alpha);
  double alphalo = alpha_min + ialpha * d_alpha;
  int iRcc = floor((Rcc - Rcc_min) / d_Rcc);
  double Rcclo = Rcc_min + iRcc * d_Rcc;
  if(ik < 0 || ik > Nk || ialpha < 0 || ialpha > Nalpha || iRcc < 0 || iRcc > NRcc)
    return 0;
  if(ik < Nk && ialpha < Nalpha && iRcc < NRcc)
  { 
    double klo_lint_iRcc = CC3_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC3_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
    double khi_lint_iRcc = CC3_array[ialpha + 1][ik][iRcc] * (klo - k + d_k) / d_k + CC3_array[ialpha + 1][ik + 1][iRcc] * (k - klo) /d_k;
    double klo_lint_iRcc_plus = CC3_array[ialpha][ik][iRcc + 1] * (klo - k + d_k) / d_k + CC3_array[ialpha][ik + 1][iRcc + 1] * (k - klo) / d_k;
    double khi_lint_iRcc_plus = CC3_array[ialpha + 1][ik][iRcc + 1] * (klo - k + d_k) / d_k + CC3_array[ialpha + 1][ik + 1][iRcc + 1] * (k - klo) /d_k;

    double alphalo_lint = klo_lint_iRcc * (alphalo - alpha + d_alpha) / d_alpha + khi_lint_iRcc * (alpha - alphalo) / d_alpha;
    double alphahi_lint = klo_lint_iRcc_plus * (alphalo - alpha + d_alpha) / d_alpha + khi_lint_iRcc_plus * (alpha - alphalo) / d_alpha;
    return alphalo_lint * (Rcclo - Rcc + d_Rcc) / d_Rcc + alphahi_lint * (Rcc - Rcclo) / d_Rcc;
  }
  else if(ik == Nk && ialpha < Nalpha && iRcc < NRcc)
  {
    double Rcclo_lint = CC3_array[ialpha][ik][iRcc] * (Rcclo - Rcc + d_Rcc) / d_Rcc + CC3_array[ialpha][ik][iRcc + 1] * (Rcc - Rcclo) / d_Rcc;
    double Rcchi_lint = CC3_array[ialpha + 1][ik][iRcc] * (Rcclo - Rcc + d_Rcc) / d_Rcc + CC3_array[ialpha + 1][ik][iRcc + 1] * (Rcc - Rcclo) /d_Rcc;
    return Rcclo_lint * (alphalo - alpha + d_alpha) / d_alpha + Rcchi_lint * (alpha - alphalo) / d_alpha;
  }
  else if(ialpha == Nalpha && ik < Nk && iRcc < NRcc)
  {
    double klo_lint = CC3_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC3_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
    double khi_lint = CC3_array[ialpha][ik][iRcc + 1] * (klo - k + d_k) / d_k + CC3_array[ialpha][ik + 1][iRcc + 1] * (k - klo) /d_k;
    return klo_lint * (Rcclo - Rcc + d_Rcc) / d_Rcc + khi_lint * (Rcc - Rcclo) / d_Rcc;
  }
  else if(iRcc == NRcc && ik < Nk && ialpha < Nalpha)
  {
    double klo_lint = CC3_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC3_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
    double khi_lint = CC3_array[ialpha + 1][ik][iRcc] * (klo - k + d_k) / d_k + CC3_array[ialpha + 1][ik + 1][iRcc] * (k - klo) /d_k;
    return klo_lint * (alphalo - alpha + d_alpha) / d_alpha + khi_lint * (alpha - alphalo) / d_alpha;
  }
  else if(iRcc == NRcc && ik == Nk && ialpha < Nalpha)
    return CC3_array[ialpha][ik][iRcc] * (alphalo - alpha + d_alpha) / d_alpha + CC3_array[ialpha + 1][ik][iRcc] * (alpha - alphalo) / d_alpha;
  else if(iRcc == NRcc && ik < Nk && ialpha == Nalpha)
    return CC3_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC3_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
  else if(iRcc < NRcc && ik == Nk && ialpha == Nalpha)
    return CC3_array[ialpha][ik][iRcc] * (Rcclo - Rcc + d_Rcc) / d_Rcc + CC3_array[ialpha][ik][iRcc + 1] * (Rcc - Rcclo) / d_Rcc;
  else if(ik == Nk && ialpha == Nalpha && iRcc == NRcc)
    return CC3_array[ialpha][ik][iRcc];
  return 0;
}

double StrongCorrFunc_reader::get_value_CC4_real(const double alpha, const double k, const double Rcc) const
{
  if(!CC4_real_array)
    return 0;
  int ik = floor((k - k_min) / d_k);
  double klo = k_min + ik * d_k;
  int ialpha = floor((alpha - alpha_min) / d_alpha);
  double alphalo = alpha_min + ialpha * d_alpha;
  int iRcc = floor((Rcc - Rcc_min) / d_Rcc);
  double Rcclo = Rcc_min + iRcc * d_Rcc;
  if(ik < 0 || ik > Nk || ialpha < 0 || ialpha > Nalpha || iRcc < 0 || iRcc > NRcc)
    return 0;
  if(ik < Nk && ialpha < Nalpha && iRcc < NRcc)
  { 
    double klo_lint_iRcc = CC4_real_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC4_real_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
    double khi_lint_iRcc = CC4_real_array[ialpha + 1][ik][iRcc] * (klo - k + d_k) / d_k + CC4_real_array[ialpha + 1][ik + 1][iRcc] * (k - klo) /d_k;
    double klo_lint_iRcc_plus = CC4_real_array[ialpha][ik][iRcc + 1] * (klo - k + d_k) / d_k + CC4_real_array[ialpha][ik + 1][iRcc + 1] * (k - klo) / d_k;
    double khi_lint_iRcc_plus = CC4_real_array[ialpha + 1][ik][iRcc + 1] * (klo - k + d_k) / d_k + CC4_real_array[ialpha + 1][ik + 1][iRcc + 1] * (k - klo) /d_k;

    double alphalo_lint = klo_lint_iRcc * (alphalo - alpha + d_alpha) / d_alpha + khi_lint_iRcc * (alpha - alphalo) / d_alpha;
    double alphahi_lint = klo_lint_iRcc_plus * (alphalo - alpha + d_alpha) / d_alpha + khi_lint_iRcc_plus * (alpha - alphalo) / d_alpha;
    return alphalo_lint * (Rcclo - Rcc + d_Rcc) / d_Rcc + alphahi_lint * (Rcc - Rcclo) / d_Rcc;
  }
  else if(ik == Nk && ialpha < Nalpha && iRcc < NRcc)
  {
    double Rcclo_lint = CC4_real_array[ialpha][ik][iRcc] * (Rcclo - Rcc + d_Rcc) / d_Rcc + CC4_real_array[ialpha][ik][iRcc + 1] * (Rcc - Rcclo) / d_Rcc;
    double Rcchi_lint = CC4_real_array[ialpha + 1][ik][iRcc] * (Rcclo - Rcc + d_Rcc) / d_Rcc + CC4_real_array[ialpha + 1][ik][iRcc + 1] * (Rcc - Rcclo) /d_Rcc;
    return Rcclo_lint * (alphalo - alpha + d_alpha) / d_alpha + Rcchi_lint * (alpha - alphalo) / d_alpha;
  }
  else if(ialpha == Nalpha && ik < Nk && iRcc < NRcc)
  {
    double klo_lint = CC4_real_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC4_real_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
    double khi_lint = CC4_real_array[ialpha][ik][iRcc + 1] * (klo - k + d_k) / d_k + CC4_real_array[ialpha][ik + 1][iRcc + 1] * (k - klo) /d_k;
    return klo_lint * (Rcclo - Rcc + d_Rcc) / d_Rcc + khi_lint * (Rcc - Rcclo) / d_Rcc;
  }
  else if(iRcc == NRcc && ik < Nk && ialpha < Nalpha)
  {
    double klo_lint = CC4_real_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC4_real_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
    double khi_lint = CC4_real_array[ialpha + 1][ik][iRcc] * (klo - k + d_k) / d_k + CC4_real_array[ialpha + 1][ik + 1][iRcc] * (k - klo) /d_k;
    return klo_lint * (alphalo - alpha + d_alpha) / d_alpha + khi_lint * (alpha - alphalo) / d_alpha;
  }
  else if(iRcc == NRcc && ik == Nk && ialpha < Nalpha)
    return CC4_real_array[ialpha][ik][iRcc] * (alphalo - alpha + d_alpha) / d_alpha + CC4_real_array[ialpha + 1][ik][iRcc] * (alpha - alphalo) / d_alpha;
  else if(iRcc == NRcc && ik < Nk && ialpha == Nalpha)
    return CC4_real_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC4_real_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
  else if(iRcc < NRcc && ik == Nk && ialpha == Nalpha)
    return CC4_real_array[ialpha][ik][iRcc] * (Rcclo - Rcc + d_Rcc) / d_Rcc + CC4_real_array[ialpha][ik][iRcc + 1] * (Rcc - Rcclo) / d_Rcc;
  else if(ik == Nk && ialpha == Nalpha && iRcc == NRcc)
    return CC4_real_array[ialpha][ik][iRcc];
  return 0;
}

double StrongCorrFunc_reader::get_value_CC4_imag(const double alpha, const double k, const double Rcc) const
{
  if(!CC4_imag_array)
    return 0;
  int ik = floor((k - k_min) / d_k);
  double klo = k_min + ik * d_k;
  int ialpha = floor((alpha - alpha_min) / d_alpha);
  double alphalo = alpha_min + ialpha * d_alpha;
  int iRcc = floor((Rcc - Rcc_min) / d_Rcc);
  double Rcclo = Rcc_min + iRcc * d_Rcc;
  if(ik < 0 || ik > Nk || ialpha < 0 || ialpha > Nalpha || iRcc < 0 || iRcc > NRcc)
    return 0;
  if(ik < Nk && ialpha < Nalpha && iRcc < NRcc)
  { 
    double klo_lint_iRcc = CC4_imag_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC4_imag_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
    double khi_lint_iRcc = CC4_imag_array[ialpha + 1][ik][iRcc] * (klo - k + d_k) / d_k + CC4_imag_array[ialpha + 1][ik + 1][iRcc] * (k - klo) /d_k;
    double klo_lint_iRcc_plus = CC4_imag_array[ialpha][ik][iRcc + 1] * (klo - k + d_k) / d_k + CC4_imag_array[ialpha][ik + 1][iRcc + 1] * (k - klo) / d_k;
    double khi_lint_iRcc_plus = CC4_imag_array[ialpha + 1][ik][iRcc + 1] * (klo - k + d_k) / d_k + CC4_imag_array[ialpha + 1][ik + 1][iRcc + 1] * (k - klo) /d_k;

    double alphalo_lint = klo_lint_iRcc * (alphalo - alpha + d_alpha) / d_alpha + khi_lint_iRcc * (alpha - alphalo) / d_alpha;
    double alphahi_lint = klo_lint_iRcc_plus * (alphalo - alpha + d_alpha) / d_alpha + khi_lint_iRcc_plus * (alpha - alphalo) / d_alpha;
    return alphalo_lint * (Rcclo - Rcc + d_Rcc) / d_Rcc + alphahi_lint * (Rcc - Rcclo) / d_Rcc;
  }
  else if(ik == Nk && ialpha < Nalpha && iRcc < NRcc)
  {
    double Rcclo_lint = CC4_imag_array[ialpha][ik][iRcc] * (Rcclo - Rcc + d_Rcc) / d_Rcc + CC4_imag_array[ialpha][ik][iRcc + 1] * (Rcc - Rcclo) / d_Rcc;
    double Rcchi_lint = CC4_imag_array[ialpha + 1][ik][iRcc] * (Rcclo - Rcc + d_Rcc) / d_Rcc + CC4_imag_array[ialpha + 1][ik][iRcc + 1] * (Rcc - Rcclo) /d_Rcc;
    return Rcclo_lint * (alphalo - alpha + d_alpha) / d_alpha + Rcchi_lint * (alpha - alphalo) / d_alpha;
  }
  else if(ialpha == Nalpha && ik < Nk && iRcc < NRcc)
  {
    double klo_lint = CC4_imag_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC4_imag_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
    double khi_lint = CC4_imag_array[ialpha][ik][iRcc + 1] * (klo - k + d_k) / d_k + CC4_imag_array[ialpha][ik + 1][iRcc + 1] * (k - klo) /d_k;
    return klo_lint * (Rcclo - Rcc + d_Rcc) / d_Rcc + khi_lint * (Rcc - Rcclo) / d_Rcc;
  }
  else if(iRcc == NRcc && ik < Nk && ialpha < Nalpha)
  {
    double klo_lint = CC4_imag_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC4_imag_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
    double khi_lint = CC4_imag_array[ialpha + 1][ik][iRcc] * (klo - k + d_k) / d_k + CC4_imag_array[ialpha + 1][ik + 1][iRcc] * (k - klo) /d_k;
    return klo_lint * (alphalo - alpha + d_alpha) / d_alpha + khi_lint * (alpha - alphalo) / d_alpha;
  }
  else if(iRcc == NRcc && ik == Nk && ialpha < Nalpha)
    return CC4_imag_array[ialpha][ik][iRcc] * (alphalo - alpha + d_alpha) / d_alpha + CC4_imag_array[ialpha + 1][ik][iRcc] * (alpha - alphalo) / d_alpha;
  else if(iRcc == NRcc && ik < Nk && ialpha == Nalpha)
    return CC4_imag_array[ialpha][ik][iRcc] * (klo - k + d_k) / d_k + CC4_imag_array[ialpha][ik + 1][iRcc] * (k - klo) / d_k;
  else if(iRcc < NRcc && ik == Nk && ialpha == Nalpha)
    return CC4_imag_array[ialpha][ik][iRcc] * (Rcclo - Rcc + d_Rcc) / d_Rcc + CC4_imag_array[ialpha][ik][iRcc + 1] * (Rcc - Rcclo) / d_Rcc;
  else if(ik == Nk && ialpha == Nalpha && iRcc == NRcc)
    return CC4_imag_array[ialpha][ik][iRcc];
  return 0;
}

void StrongCorrFunc_reader::InitTable(const char* filename)
{ 
  ifstream infile;
  infile.open(filename, ios::in | ios::binary);
  cout << "Initializing StrongCorrFunc_reader table in memory using file: " << filename << endl;
  infile.read((char*)&alpha_min, sizeof(alpha_min));  cout << "  alpha_min = " << alpha_min << endl;
  infile.read((char*)&alpha_max, sizeof(alpha_max));  cout << "  alpha_max = " << alpha_max << endl;
  infile.read((char*)&Nalpha,    sizeof(Nalpha));     cout << "  Nalpha = " << Nalpha << endl;
  infile.read((char*)&k_min,     sizeof(k_min));      cout << "  k_min = " << k_min << endl;
  infile.read((char*)&k_max,     sizeof(k_max));      cout << "  k_max = " << k_max << endl;
  infile.read((char*)&Nk,        sizeof(Nk));         cout << "  Nk = " << Nk << endl;
  infile.read((char*)&Rcc_min,   sizeof(Rcc_min));    cout << "  Rcc_min = " << Rcc_min << endl;
  infile.read((char*)&Rcc_max,   sizeof(Rcc_max));    cout << "  Rcc_max = " << Rcc_max << endl;
  infile.read((char*)&NRcc,      sizeof(NRcc));       cout << "  NRcc = " << NRcc << endl;
  d_alpha = (alpha_max - alpha_min) * 1.0 / Nalpha;
  d_k = (k_max - k_min) * 1.0 / Nk;
  d_Rcc= (Rcc_max - Rcc_min) * 1.0 / NRcc;
  cout << "    d_alpha: " << d_alpha  << " d_k: " << d_k << " d_Rcc: " << d_Rcc << endl;
  double value1 = 0., value2 = 0., value3 = 0., value4_real = 0., value4_imag = 0.;
//  Nalpha-=10;
//  cout << "Nalpha: " << Nalpha << endl;
  CC1_array = new double**[Nalpha + 1];
  CC2_array = new double**[Nalpha + 1];
  CC3_array = new double**[Nalpha + 1];
  CC4_real_array = new double**[Nalpha + 1];
  CC4_imag_array = new double**[Nalpha + 1];
  for(int ialpha = 0; ialpha < (Nalpha + 1); ialpha++)
  {
    cout << ialpha << endl;
    CC1_array[ialpha] = new double*[Nk + 1];
    CC2_array[ialpha] = new double*[Nk + 1];
    CC3_array[ialpha] = new double*[Nk + 1];
    CC4_real_array[ialpha] = new double*[Nk + 1];
    CC4_imag_array[ialpha] = new double*[Nk + 1];
    for(int ik = 0; ik < (Nk + 1); ik++)
    {
      CC1_array[ialpha][ik] = new double[NRcc + 1];
      CC2_array[ialpha][ik] = new double[NRcc + 1];
      CC3_array[ialpha][ik] = new double[NRcc + 1];
      CC4_real_array[ialpha][ik] = new double[NRcc + 1];
      CC4_imag_array[ialpha][ik] = new double[NRcc + 1];
    }
  }

  for(int ialpha = 0; ialpha < (Nalpha + 1); ialpha++)
  {
    for(int ik = 0; ik < (Nk + 1); ik++)
    {
      for(int iRcc = 0; iRcc < (NRcc + 1); iRcc++)
      {
        infile.read((char*)&value1, sizeof(value1));
        infile.read((char*)&value2, sizeof(value2));
        infile.read((char*)&value3, sizeof(value3));
        infile.read((char*)&value4_real, sizeof(value4_real));
        infile.read((char*)&value4_imag, sizeof(value4_imag));
        CC1_array[ialpha][ik][iRcc] = value1;
        CC2_array[ialpha][ik][iRcc] = value2;
        CC3_array[ialpha][ik][iRcc] = value3;
        CC4_real_array[ialpha][ik][iRcc] = value4_real;
        CC4_imag_array[ialpha][ik][iRcc] = value4_imag;
      }
    }
  }
  cout << "CorrFunc_table loaded." << endl;
  infile.close();
}
