#include "../../Torus/src/utils/Units.h"
using namespace Units;
namespace py = pybind11;
void init_units(py::module_ &torus) {
	torus.attr("angle_unit") = angle_unit;
	torus.attr("length_unit") = length_unit;
	torus.attr("time_unit") = time_unit;
	torus.attr("mass_unit") = mass_unit;
		// angle
	torus.attr("radian") = radian;
	torus.attr("kpc") = kpc;
	torus.attr("Myr") = Myr;
	torus.attr("Msun") = Msun;
	torus.attr("rad") = rad;
	torus.attr("degree") = degree;
	torus.attr("arcmin") = arcmin;
	torus.attr("angsec") = angsec;
		// length
	torus.attr("mas") = mas;
	torus.attr("anghr") = anghr;
	torus.attr("angmin") = angmin;
	torus.attr("angsec") = angsec;
	torus.attr("cm") = cm;
	torus.attr("meter") = meter;
	torus.attr("km") = km;
	torus.attr("ly") = ly;
	torus.attr("pc") = pc;
	torus.attr("Mpc") = Mpc;
	torus.attr("AU") = AU;
		// time
	torus.attr("sec") = sec;
	torus.attr("hour") = hour;
	torus.attr("day") = day;
	torus.attr("yr") = yr;
	torus.attr("hyr") = hyr;
	torus.attr("century") = century;
	torus.attr("Gyr") = Gyr;
		// mass
	torus.attr("gram") = gram;
	torus.attr("kg") = kg;
		// velocity
	torus.attr("kpcMyr") = kpcMyr;
	torus.attr("kpcGyr") = kpcGyr;
	torus.attr("AUyr") = AUyr;
	torus.attr("kms") = kms;
	torus.attr("c_light") = c_light;
		// angle velocity
	torus.attr("radMyr") = radMyr;
	torus.attr("kmskpc") = kmskpc;
	torus.attr("masyr") = masyr;
	torus.attr("ashyr") = ashyr;
	torus.attr("secyr") = secyr;
	torus.attr("asyr") = asyr;
		// area
	torus.attr("pc2") = pc2;
	torus.attr("kpc2") = kpc2;
	torus.attr("cm2") = cm2;
		// volume
	torus.attr("pc3") = pc3;
	torus.attr("kpc3") = kpc3;
	torus.attr("cm3") = cm3;
		// constant of gravity
	torus.attr("G") = G;
	torus.attr("Grav") = Grav;
	torus.attr("fPiG") = fPiG;
		// inverse quantities 
		// inverse angle
	torus.attr("rad_i") = rad_i;
	torus.attr("radian_i") = radian_i;
	torus.attr("degree_i") = degree_i;
	torus.attr("arcmin_i") = arcmin_i;
	torus.attr("arcsec_i") = arcsec_i;
	torus.attr("mas_i") = mas_i;
	torus.attr("anghr_i") = anghr_i;
	torus.attr("angmin_i") = angmin_i;
	torus.attr("angsec_i") = angsec_i;
		// inverse length
	torus.attr("cm_i") = cm_i;
	torus.attr("meter_i") = meter_i;
	torus.attr("km_i") = km_i;
	torus.attr("AU_i") = AU_i;
	torus.attr("ly_i") = ly_i;
	torus.attr("pc_i") = pc_i;
	torus.attr("kpc_i") = kpc_i;
	torus.attr("Mpc_i") = Mpc_i;
		// inverse time
	torus.attr("sec_i") = sec_i;
	torus.attr("hour_i") = hour_i;
	torus.attr("day_i") = day_i;
	torus.attr("yr_i") = yr_i;
	torus.attr("hyr_i") = hyr_i;
	torus.attr("century_i") = century_i;
	torus.attr("Myr_i") = Myr_i;
	torus.attr("Gyr_i") = Gyr_i;
		// inverse mass
	torus.attr("gram_i") = gram_i;
	torus.attr("kg_i") = kg_i;
	torus.attr("Msun_i") = Msun_i;
		// inverse velocity
	torus.attr("kpcMyr_i") = kpcMyr_i;
	torus.attr("kpcGyr_i") = kpcGyr_i;
	torus.attr("AUyr_i") = AUyr_i;
	torus.attr("kms_i") = kms_i;
	torus.attr("c_light_i") = c_light_i;
		// inverse angle velocity
	torus.attr("radMyr_i") = radMyr_i;
	torus.attr("kmskpc_i") = kmskpc_i;
	torus.attr("masyr_i") = masyr_i;
	torus.attr("ashyr_i") = ashyr_i;
	torus.attr("secyr_i") = secyr_i;
	torus.attr("asyr_i") = asyr_i;
		// inverse area
	torus.attr("pc2_i") = pc2_i;
	torus.attr("kpc2_i") = kpc2_i;
	torus.attr("cm2_i") = cm2_i;
		// inverse volume
	torus.attr("pc3_i") = pc3_i;
	torus.attr("kpc3_i") = kpc3_i;
	torus.attr("cm3_i") = cm3_i;
		// inverse of constant of gravity
	torus.attr("G_i") = G_i;
	torus.attr("Grav_i") = Grav_i;
	torus.attr("fPiG_i") = fPiG_i;
}
