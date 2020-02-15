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
	double a, da, e, de, i, di, Om, dOm, w, dw, Tp, dTp, P, dP, Kmag, R_, dR_, q_ang, dq_ang;
	int ref;

	// First output the format-dependent header
	switch(format) {
		case(ascii):
			fprintf(ofile, "%-4s  %6s %6s  %6s %6s  %6s %4s  %6s %5s  %6s %5s  %8s %7s  %7s %6s  %5s   %8s %8s %6s %6s  %6s\n", "id", "a", "da", "e", "de", "i", "di", "Om", "dOm", "w", "dw", "Tp", "dTp", "P", "dP", "Kmag", "q", "dq", "v", "dv", "vavg");
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
				"\t\t<th>v&#773;</th>\n"
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
				"! &Delta;v\n"
				"! v&#773;\n");
			break;
	}

	for(int iline = 0; (nread=getline(&line, &bufsize, ifile)) != -1; iline++) {
		if(nread > 0 && line[0] == '#') continue;
		nscan = sscanf(line, "%ms %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf",
				&name, &a, &da, &e, &de, &i, &di, &Om, &dOm, &w, &dw, &Tp, &dTp, &P, &dP, &Kmag, &R_, &dR_, &ref, &q_ang, &dq_ang);
		if(nscan != 21) {
			fprintf(stderr, "Error parsing line %d starting with: %.15s\n", iline+1, line);
			continue;
		}

		// 1st order analytical error bars, ignoring correlations, but using r_ang
		if(q_ang < 0) q_ang = a*(1-e);
		double q       = R*q_ang;
		double dq      = q*sqrt(dR*dR/(R*R)+dq_ang*dq_ang/(q_ang*q_ang));
		double v       = 2*pi*a*R/P*sqrt(2*a/q_ang-1);
		double dv      = 2*pi*a*R/P*sqrt( (2*a/q_ang-1)*(da*da/(a*a)+dR*dR/(R*R)+dP*dP/(P*P)) + dq_ang*dq_ang/(4*q_ang*q_ang) + da*da/(a*a)*(sqrt(2*a/q_ang-1)+0.25) );

		// Get mean velocity too. This is based on
		// https://www.mathsisfun.com/geometry/ellipse-perimeter.html approx #3 by Ramanujan
		// It's hard to ge an accurate uncertainty on this. It should be smaller than dv, though.
		double b = a*sqrt(1-e*e), h = (a-b)*(a-b)/((a+b)*(a+b));
		double circum  = pi*(a+b)*(1+3*h/(10+sqrt(4-3*h))) * R;
		double v_avg   = circum/P * pc*arcsec/yr/c;

		// Convert to distance in AU and speed in beta
		q *= pc*arcsec/AU;   dq *= pc*arcsec/AU;
		v *= pc*arcsec/yr/c; dv *= pc*arcsec/yr/c;

		// Output table row
		switch(format) {
			case(ascii):
				fprintf(ofile, "%-4s  %6.4f %6.4f  %6.4f %6.4f  %6.2f %4.2f  %6.2f %5.2f  %6.2f %5.2f  %8.3f %7.3f  %7.1f %6.1f  %5.2f   %8.1f %8.1f %6.2f %6.2f  %6.2f\n", name, a, da, e, de, i, di, Om, dOm, w, dw, Tp, dTp, P, dP, Kmag,  q, dq, v*100, dv*100, v_avg*100);
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
					"\t\t<td>%6.2f</td>\n"
					"\t</tr>\n", name, a, da, e, de, i, di, Om, dOm, w, dw, Tp, dTp, P, dP, Kmag,  q, dq, v*100, dv*100, v_avg*100);
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
					"| %6.2f\n"
					"| %6.2f\n", name, a, da, e, de, i, di, Om, dOm, w, dw, Tp, dTp, P, dP, Kmag,  q, dq, v*100, dv*100, v_avg*100);
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
