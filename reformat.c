#define _POSIX_C_SOURCE 201910L

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>

#define pi 3.14159265359
#define yr 31556925.216
#define AU 149597870700e0
#define pc 3.08567758149137e16
#define c  299792458e0
#define arcsec 4.84813681109536e-06

enum { ascii, html, wiki };

void help() {
	fprintf(stderr, "Usage: reformat [-f format] [params.txt] [ofile], where format can be ascii, html or wiki\n");
	exit(1);
}

int main(int argc, char ** argv) {
	char * ifname = "/dev/stdin";
	char * ofname = "/dev/stdout";
	int format = ascii;
	for(int i = 1, narg = 0; i < argc; i++) {
		if(!strcmp(argv[i], "-f")) {
			if(++i >= argc) help();
			if     (!strcmp(argv[i], "ascii")) format = ascii;
			else if(!strcmp(argv[i], "html"))  format = html;
			else if(!strcmp(argv[i], "wiki"))  format = wiki;
			else help();
		}
		else if(strlen(argv[i]) > 0 && argv[i][0] == '-') help();
		else {
			if     (narg == 0) ifname = argv[i];
			else if(narg == 1) ofname = argv[i];
			else help();
			narg++;
		}
	}
	FILE * ifile = fopen(ifname, "r");
	if(!ifile) { fprintf(stderr, "Error opening '%s' for reading!\n", ifname); exit(1); }
	FILE * ofile = fopen(ofname, "w");
	if(!ofile) { fprintf(stderr, "Error opening '%s' for writing!\n", ofname); exit(1); }
	char * line = NULL, *name;
	ssize_t nread, nscan; size_t bufsize;
	int nsamp = 10000, si;

	double R = 8179, dR = 13;
	double a, da, e, de, i, di, Om, dOm, w, dw, Tp, dTp, P, dP, Kmag, R_, dR_, dr_ang;
	int ref;

	// First output the format-dependent header
	switch(format) {
		case(ascii):
			fprintf(ofile, "%-4s  %6s %6s  %6s %6s  %6s %4s  %6s %5s  %6s %5s  %8s %7s  %7s %6s  %5s   %8s %8s %6s %6s\n", "id", "a", "da", "e", "de", "i", "di", "Om", "dOm", "w", "dw", "Tp", "dTp", "P", "dP", "Kmag", "q", "dq", "v", "dv");
			break;
		case(html):
			fprintf(ofile, "<table>\n\t<tr>\n"
				"\t\t<th>id</th><th>a (\")</th>\n"
				"\t\t<th>&Delta;a</th>\n"
				"\t\t<th>e</th>\n"
				"\t\t<th>&Delta;e</th>\n"
				"\t\t<th>i (°)</th>\n"
				"\t\t<th>&Delta;i</th>\n"
				"\t\t<th>&Omega; (°)</th>\n"
				"\t\t<th>&Delta;&Omega;</th>\n"
				"\t\t<th>&omega; (°)</th>\n"
				"\t\t<th>&Delta;&omega;</th>\n"
				"\t\t<th>Tp (yr)</th>\n"
				"\t\t<th>&Delta;Tp</th>\n"
				"\t\t<th>P (yr)</th>\n"
				"\t\t<th>&Delta;P</th>\n"
				"\t\t<th>Kmag</th>\n"
				"\t\t<th>q (AU)</th>\n"
				"\t\t<th>&Delta;q</th>\n"
				"\t\t<th>v (%%c)</th>\n"
				"\t\t<th>&Delta;v</th>\n"
				"\t</tr>\n");
			break;
		case(wiki):
			fprintf(ofile, "{|class=\"wikitable sortable\"  style=\"text-align:right\"\n"
				"! id\n"
				"! a\n"
				"! &Delta;a\n"
				"! e\n"
				"! &Delta;e\n"
				"! i (°)\n"
				"! &Delta;i\n"
				"! &Omega; (°)\n"
				"! &Delta;&Omega;\n"
				"! &omega; (°)\n"
				"! &Delta;&omega;\n"
				"! Tp (yr)\n"
				"! &Delta;Tp\n"
				"! P (yr)\n"
				"! &Delta;P\n"
				"! Kmag\n"
				"! q (AU)\n"
				"! &Delta;q\n"
				"! v (%%c)\n"
				"! &Delta;v\n");
			break;
	}

	for(int iline = 0; (nread=getline(&line, &bufsize, ifile)) != -1; iline++) {
		if(nread > 0 && line[0] == '#') continue;
		nscan = sscanf(line, "%ms %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf",
				&name, &a, &da, &e, &de, &i, &di, &Om, &dOm, &w, &dw, &Tp, &dTp, &P, &dP, &Kmag, &R_, &dR_, &ref, &dr_ang);
		if(nscan != 20) {
			fprintf(stderr, "Error parsing line %d starting with: %.15s\n", iline+1, line);
			continue;
		}

		// 1st order analytical error bars, ignoring correlations, but using r_ang
		double r_ang   = a*(1-e);
		double r       = R*r_ang;
		double dr      = r*sqrt(dR*dR/(R*R)+dr_ang*dr_ang/(r_ang*r_ang));
		double v       = 2*pi*a*R/P*sqrt(2*a/r_ang-1);
		double dv      = 2*pi*a*R/P*sqrt( (2*a/r_ang-1)*(da*da/(a*a)+dR*dR/(R*R)+dP*dP/(P*P)) + dr_ang*dr_ang/(4*r_ang*r_ang) + da*da/(a*a)*(sqrt(2*a/r_ang-1)+0.25) );

		// Convert to distance in AU and speed in beta
		r *= pc*arcsec/AU;   dr *= pc*arcsec/AU;
		v *= pc*arcsec/yr/c; dv *= pc*arcsec/yr/c;

		// Output table row
		switch(format) {
			case(ascii):
				fprintf(ofile, "%-4s  %6.4f %6.4f  %6.4f %6.4f  %6.2f %4.2f  %6.2f %5.2f  %6.2f %5.2f  %8.3f %7.3f  %7.1f %6.1f  %5.2f   %8.1f %8.1f %6.2f %6.2f\n", name, a, da, e, de, i, di, Om, dOm, w, dw, Tp, dTp, P, dP, Kmag,  r, dr, v*100, dv*100);
				break;
			case(html):
				fprintf(ofile, "\t<tr>\n"
					"\t\t<td>%s</td>\n"
					"\t\t<td>%6.4f</td>\n"
					"\t\t<td>%6.4f</td>\n"
					"\t\t<td>%6.4f</td>\n"
					"\t\t<td>%6.4f</td>\n"
					"\t\t<td>%6.2f</td>\n"
					"\t\t<td>%4.2f</td>\n"
					"\t\t<td>%6.2f</td>\n"
					"\t\t<td>%5.2f</td>\n"
					"\t\t<td>%6.2f</td>\n"
					"\t\t<td>%5.2f</td>\n"
					"\t\t<td>%8.3f</td>\n"
					"\t\t<td>%7.3f</td>\n"
					"\t\t<td>%7.1f</td>\n"
					"\t\t<td>%6.1f</td>\n"
					"\t\t<td>%5.2f</td>\n"
					"\t\t<td>%8.1f</td>\n"
					"\t\t<td>%8.1f</td>\n"
					"\t\t<td>%6.2f</td>\n"
					"\t\t<td>%6.2f</td>\n"
					"\t</tr>\n", name, a, da, e, de, i, di, Om, dOm, w, dw, Tp, dTp, P, dP, Kmag,  r, dr, v*100, dv*100);
				break;
			case(wiki):
				fprintf(ofile, "|-\n"
					"| %s\n"
					"| %6.4f\n"
					"| %6.4f\n"
					"| %6.4f\n"
					"| %6.4f\n"
					"| %6.2f\n"
					"| %4.2f\n"
					"| %6.2f\n"
					"| %5.2f\n"
					"| %6.2f\n"
					"| %5.2f\n"
					"| %8.3f\n"
					"| %7.3f\n"
					"| %7.1f\n"
					"| %6.1f\n"
					"| %5.2f\n"
					"| %8.1f\n"
					"| %8.1f\n"
					"| %6.2f\n"
					"| %6.2f\n", name, a, da, e, de, i, di, Om, dOm, w, dw, Tp, dTp, P, dP, Kmag,  r, dr, v*100, dv*100);
				break;
		}
		free(name);
	}
	// Output footer
	switch(format) {
		case(ascii): break;
		case(html): fprintf(ofile, "</table>\n"); break;
		case(wiki): fprintf(ofile, "|}\n"); break;
	}

	free(line);
	fclose(ifile);
	fclose(ofile);
	return 0;
}

#if 0
// Old stuff
//
		// Compute mean and error bars, assuming independent errors, via random sampling
		double r = 0, dr = 0, v = 0, dv = 0;
		for(si = 0; si < nsamp; si++) {
			double R_ = R + 0*rand_gauss()*dR;
			double a_ = a + rand_gauss()*da;
			double e_ = fmin(e + rand_gauss()*de, 1-1e-3);
			double P_ = P + rand_gauss()*dP;
			double r_ = R_*a_*(1-e_);
			double v_ = 2*pi*R_*a_/P_*sqrt(1-e_*e_)/(1-e_);
			r += r_; dr += r_*r_;
			v += v_; dv += v_*v_;
		}
		r /= nsamp; dr = sqrt(dr/nsamp-r*r); // pc*arcsec
		v /= nsamp; dv = sqrt(dv/nsamp-v*v); // pc*arcsec/yr

// Random number stuff
static uint64_t       state      = 0x4d595df4d0f33173;   // Or something seed-dependent
static uint64_t const multiplier = 6364136223846793005u;
static uint64_t const increment  = 1442695040888963407u; // Or an arbitrary odd constant
static uint32_t rotr32(uint32_t x, unsigned r) { return x >> r | x << (-r & 31); }
uint32_t pcg32(void)
{
	uint64_t x = state;
	unsigned count = (unsigned)(x >> 59);      // 59 = 64 - 5
	state = x * multiplier + increment;
	x ^= x >> 18;                              // 18 = (64 - 27)/2
	return rotr32((uint32_t)(x >> 27), count); // 27 = 32 - 5
}
void pcg32_init(uint64_t seed)
{
	state = seed + increment;
	(void)pcg32();
}

double rand_uniform() {
	// map directly onto float in range [1,2), then subtract 1
	uint64_t rint = pcg32();
	uint64_t work = (0x3fflu<<52)|(rint<<20);
	return (*(double*)&work)-1;
}

double rand_gauss() {
	double epsilon = 2e-314;
	double u1, u2;
	do {
		u1 = rand_uniform();
		u2 = rand_uniform();
	}
	while (u1 <= epsilon);
	return sqrt(-2.0*log(u1))*cos(2*pi*u2);
}

#endif
