//
// Single parenthood and adolescent sex-related behaviors
//

// Dufur, Hoffmann, and Erickson

clear
set more off
set emptycells drop
set matsize 1000
cd "/Volumes/Add Health/dta"

//////////////
//          //
//   DATA   //
//          //
//////////////

// wave I
use aid bio_sex h1gh26 h1gh24 h1gh25 pa1 pa2 pa4 pa6_* pa38 pb2 pb8 pa8b pa10 pa11 pa12 pa17 pb13 pa38-pa44 pa55 pc1 pc6b h1hr* h1rm7 h1rf7 h1rm1 h1rf1 h1gi6* h1gi4 h1rp* h1se1-h1se3 h1mo* h1kq* h1id1f h1id1l-h1id1q h1co* h1bc* h1ir1 h1ir2 h1ir3 h1nr3 h1nr5 h1nr6 h1pa* h1da1 h1da3 h1da4 h1da5 h1da6 h1da8 h1gh1 h1gh1a h1gh23a h1gh23g h1gh23j h1gh28 h1gh30* h1gh31* h1gh32-h1gh44 h1gh46 h1gh48-h1gh53 h1gh44 h1gh46 h1gh57-h1gh60 h1hs1 h1hs3 h1hs5 h1hs9 using "allwave1", clear

// created wave I variables
merge 1:1 aid using "vars1", keepusing (aid agew1 famst5) nogen

// inschool	// using a 1:m match because of the missing aid variables in the inschool data--otherwise there's an error
preserve
	use aid sschlcde s12 s18 using "inschool", clear
	ren sschlcde scid
	tempfile inschool
	save "`inschool'", replace
restore
merge 1:m aid using "`inschool'", nogen
drop if aid==""

keep if bio_sex!=.						// keep only Wave 1 respondenrs

///////////////////////////
//                       //
//   DATA MANIPULATION   //
//                       //
///////////////////////////

//                //
// SINGLE PARENTS //
//                //

// responding parent's sex
recode pa1 pb2 (1=0 "Male")(2=1 "Female")(6/7=.), gen(r_fem p_fem) label(sex)

// parent marital status from parent questionnaire
recode pa10 (1=1 Single)(2=2 Married)(3=3 Widowed)(4=4 Divorced)(5=5 Separated)(6=.), gen(p_mar)

// single parent status from parent questionnaire
gen singmom = .
replace singmom = 1 if p_mar!=2 & r_fem==1
replace singmom = 0 if p_mar!=2 & r_fem==0

// unambiguous single parenthood (data from household roster and parent report)
gen spmom = .
replace spmom = 1 if r_fem==1 & (p_mar!=2 & famst5==3)
replace spmom = 0 if r_fem==0 & (p_mar!=2 & famst5==4)

	// this variable only includes parents ...
		// who reported their gender
		// who self-reported their marital status as single
		// whose marital status was determined to be single from the household roster


//          //
// CONTROLS //
//          //

// child sex
recode bio_sex (1=0)(2=1)(6/9=.), into(female)

// child age
gen age=trunc(agew1)
replace age=15 if age==14			// one case that was just below 15 slipped into the questions about sex

// child race
gen race=h1gi6a								// h1gi6a = white
replace race=2 if h1gi6b==1					// h1gi6b = black
replace race=3 if h1gi6d==1					// h1gi6d = asian
replace race=5 if h1gi6c==1					// h1gi6c = native american
replace race=4 if h1gi4==1					// h1gi4  = hispanic
recode race (0=.) (6/9=.)

// PARENTS' EDUCATION by sex and data source
// 1)nohs, 2)high school, 3)some college, 4)college
recode pa12  (1/2 10=1) (3/5=2) (6/7=3) (8/9=4) (else=0), into (paed)	// paed=parent self-report
recode pb8   (1/2 10=1) (3/5=2) (6/7=3) (8/9=4) (else=0), into (sped)	// sped=parent report on spouse
recode h1rm1 (1/2 10=1) (3/5=2) (6/7=3) (8/9=4) (else=0), into (moed1)	// moed1=child in home report on mother
recode h1rf1 (1/2 10=1) (3/5=2) (6/7=3) (8/9=4) (else=0), into (faed1)	// faed1=child in home report on father
recode s12   (1/2 10=1) (3/5=2) (6/7=3) (8/9=4) (else=0), into (moed2)	// moed2=child in school report on mother
recode s18   (1/2 10=1) (3/5=2) (6/7=3) (8/9=4) (else=0), into (faed2)	// faed2=child in school report on father
	// HIGHEST PARENT EDUCATION
	// establish the highest education within each of the reports
	gen hpasped=paed
	replace hpasped=sped if (hpasped==0 | sped>paed)					// hpasped=highest parent/spouse education from inhome parent
	gen hmofaed1=moed1
	replace hmofaed1=faed1 if (hmofaed1==0 | faed1>moed1)				// hmofaed1=highest mother/father education from inhome child
	gen hmofaed2=moed2
	replace hmofaed2=faed2 if (hmofaed2==0 | faed2>moed2)				// hmofaed2=highest mother/father education from inschool child
		// for overall highest parent educ, use parent reports first, then child in-home reports, then child in-school report
		gen hpaed=hpasped
		replace hpaed=hmofaed1 if hpaed==0
		replace hpaed=hmofaed2 if hpaed==0
		recode hpaed (0=.)

// parent income
recode pa55 (9996=.), gen(income)
gen inclog = ln(1+income)

// PARENT EMPLOYMENT STATUS
recode pa17 (7=0)(6 8 9=.), gen(full)



//////////////////////////
//     sex outcomes     //
//////////////////////////

// higher values represent perceptions of sex as risky or bad
recode h1rp* (6/9=.)

// higher values represent more surity
recode h1se1-h1se3 (5=1)(4=2)(3=3)(2=4)(1=5)(6=.a)(96/99=.)

// higher values represent more social inhibitors
recode h1mo1 (6/9=.)
recode h1mo2 h1mo4 (5=1)(4=2)(3=3)(2=4)(1=5)(else=.)
	// h1mo4 was only asked of respondents who had a resident mother

// values of 1 represent a correct response
recode h1kq1a h1kq3a h1kq4a h1kq5a h1kq6a h1kq7a h1kq9a h1kq10a (2=1)(1=0)(else=.)
recode h1kq2a h1kq8a (1=1)(2=0)(else=.)

// values of 1 would do activity
recode h1id1f h1id1l-h1id1q (1=1)(2=0)(else=.)


recode h1co1 h1co3 h1co6 h1co13 h1co15 h1co16* (7=0)(6 8 9=.)
egen sti = rowmax(h1co16*)

// higher values represent high motivations to use birth control
recode h1bc* (6/9=.)

//
recode h1nr3 h1nr5 (6/9=.)
recode h1nr6 (996/999=.)

// higher values indicate more approval
recode h1pa1-h1pa6 (6/9=.)					// no cases for single dads? makes measures irrelevant
	
recode h1pa7(6/9=.)


keep if spmom!=.				// keep only single parents

///////////////////////////////////////////////////////
// checking sexual behavior of those younger than 15 //
//
fre spmom if age<15 & !missing(female,age,race,hpaed,inclog,full)
logistic h1co1 spmom if age<15 & !missing(female,age,race,hpaed,inclog,full)   // or for spmom: 5.95
logistic h1co1 spmom female age ib1.race ib1.hpaed inclog full if age<15	   // or for spmom: 4.16
///////////////////////////////////////////////////////

drop if h1rp1==7				// drop those under 15 years old (were not asked questions about sex)

save "/Volumes/Add Health/created/FamStr SEX prepped", replace



/////////////////////////////////
//                             //
//     MULTIPLE IMPUTATION     //
//                             //
/////////////////////////////////
cd "/Volumes/Add Health/created"
use "FamStr SEX prepped", clear

mi set flong 
mi register imputed ///
	race hpaed inclog 																			/// controls
	h1rp1 h1rp2 h1se1-h1se3 h1mo1-h1mo3 h1mo5-h1mo14 h1bc1-h1bc8								/// continuous sex outcomes
	h1kq*a h1id1f h1id1l-h1id1q h1co1 h1co3 h1co6 h1co13 h1co15 sti h1pa7 						//  dichotomous sex outcomes
mi fvset base 1 race hpaed

local nscore h1rp1 h1rp2 h1se1 h1se2 h1se3 h1mo1 h1mo2 h1mo3 h1mo5 h1mo6 h1mo7 h1mo8 h1mo9 h1mo10 h1mo11 h1mo12 h1mo13 h1mo14 h1bc1 h1bc2 h1bc3 h1bc4 h1bc5 h1bc6 h1bc7 h1bc8
local n : word count `nscore'
nscore `nscore', gen(n)
mi register imputed n1-n`n'

mi impute chained ///
	(regress) inclog n1-n`n' ///
	(logit)   h1kq*a h1id1f h1id1l-h1id1q h1co1 h1co3 h1co6 h1co13 h1co15 sti h1pa7  ///
	(mlogit)  race hpaed ///
	 = female age ///
	 , add(1) burnin(300) rseed(3699) savetrace("FamStr SEX impute trace", replace) augment chaindots 
 
invnscore `nscore'

save "FamStr SEX impute diagnostic", replace

// diagnostics	 
	// check range
		mi convert wide, clear
		foreach v of varlist race hpaed inclog h1rp1 h1rp2 h1se1-h1se3 h1mo1-h1mo3 h1mo5-h1mo14 h1bc1-h1bc8 h1kq*a h1id1f h1id1l-h1id1q h1co1 h1co3 h1co6 h1co13 h1co15 sti h1pa7 {
			su `v' _1_`v'
			//fre `v' _1_`v'
		}
	// open trace data
	use "FamStr SEX impute trace", clear
	tsset iter

		foreach v of varlist *_mean {
			tsline `v', title(Trace plot of {it:`v'}) name(`v', replace)								// trace plots
			ac `v', lag(300) title(Autocorrelation plot of {it:`v'}) name(`v'_, replace)				// autocorrelation plots
		}
		drop h1co13_sd
		foreach v of varlist *_sd {
			tsline `v', title(Trace plot of {it:`v'}) name(`v', replace)								// trace plots
			ac `v', lag(300) title(Autocorrelation plot of {it:`v'}) name(`v'_, replace)				// autocorrelation plots
		}
			// diagnostics look good at 100 burnins


// RUN ACTUAL IMPUTATIONS //	SEX OUTCOMES
cd "/Volumes/Add Health/created"
use "FamStr SEX prepped", clear

set matsize 1000
mi set flong 
mi register imputed ///
	race hpaed inclog full																		/// controls
	h1rp1 h1rp2 h1se1-h1se3 h1mo1-h1mo3 h1mo5-h1mo14 h1bc1-h1bc8								/// continuous sex outcomes
	h1kq*a h1id1f h1id1l-h1id1q h1co1 h1co3 h1co6 h1co13 h1co15 sti h1pa7 						 // dichotomous sex outcomes
mi fvset base 1 race hpaed

local nscore h1rp1 h1rp2 h1se1 h1se2 h1se3 h1mo1 h1mo2 h1mo3 h1mo5 h1mo6 h1mo7 h1mo8 h1mo9 h1mo10 h1mo11 h1mo12 h1mo13 h1mo14 h1bc1 h1bc2 h1bc3 h1bc4 h1bc5 h1bc6 h1bc7 h1bc8
local n : word count `nscore'
quietly nscore `nscore', gen(n)
mi register imputed n1-n`n'

//capture log c
//log using "C:\Users\lde\Documents\DATA\Family Structure\Imputations.log", replace
mi impute chained ///
	(regress) inclog n1-n`n' ///
	(logit)   full h1kq*a h1id1f h1id1l-h1id1q h1co1 h1co3 h1co6 h1co13 h1co15 sti h1pa7  ///
	(mlogit)  race hpaed ///
	 = female age ///
	 , add(20) burnin(100) rseed(3699) augment chaindots
//log c

invnscore `nscore'
save "/Volumes/Add Health/created/FamStr SEX impute", replace



/////////////////////////////////
//                             //
//      open imputed data      //
//                             //
/////////////////////////////////
cd "/Users/lde/Desktop/Computer Backup/Users/lde/Documents/DATA/Family Structure/Sex"
use "/Volumes/Add Health/created/FamStr SEX impute", clear
//set more off

capture drop nohs hsgrad somecol colgrad white black asian natam hisp

// create dummies
tab race, gen(race_)
	rename (race_*)(white black asian natam hisp)
quietly tab hpaed, gen(rechpaed)
	rename (rechpaed*)(nohs hsgrad somecol colgrad)

capture gen x=uniform()											// generate random variable to use as dependent in regression for descriptive table


//
// checking if data reduction is possible
//
factor h1mo1-h1mo8 if _mi_m==1, factor(1) blanks(.4)
factor h1mo1-h1mo8 if _mi_m==1
rotate, blanks(.4)
alpha h1mo1 h1mo3 h1mo5-h1mo8, gen(sex_full)	// alpha = .72
alpha h1mo5 h1mo6, gen(sex_phy)					// alpha = .78
alpha h1mo1 h1mo7 h1mo8, gen(sex_soc)			// alpha = .68
alpha h1mo2-h1mo3, gen(sex_resp)				// alpha = .67

factor h1mo9-h1mo14 if _mi_m==1, blanks(.4)
alpha h1mo9-h1mo14, gen(prego)					// alpha = .72

factor h1bc* if _mi_m==1
alpha h1bc*										// alpha = .78
alpha h1bc1-h1bc5, gen(birthc1)		// alpha = .81

factor h1se1-h1se3 if _mi_m==1
alpha h1se1-h1se3, gen(birthc2)					// alpha = .66

factor h1rp*
rotate, blanks(.2)
alpha h1rp1 h1rp2 

// knowledge quiz
egen know = rowtotal(h1kq*a)



////////////////////
//
//    supplemental exploration of birth control variables
//
////////////////////
capture log close
log using "/Users/lde/Desktop/Computer Backup/Users/lde/Documents/Data/Family Structure/Birth control crosstabs.log", replace
foreach x in h1co13 h1co15 {
	foreach y in birthc1 h1bc7 h1bc8 birthc2 h1co3 h1co6 {	
		tab `y' `x' if _mi_m==1 & female==1, col nofreq
	}
}
log c




/////////////////////////////////////////
//                                     //
//     BIVARIATE DESCRIPTIVE TABLE     //
//                                     //
/////////////////////////////////////////

mat drop _all

local cont	female age white black asian natam hisp nohs hsgrad somecol colgrad inclog full
local out	h1id1f h1id1l h1id1m h1id1n h1id1o h1id1p h1id1q know sex_phy sex_soc sex_resp h1co1 birthc1 h1bc7 h1bc8 birthc2 h1co13 h1co15 h1co3 h1co6 h1rp2 prego sti

local var `out' `cont'

/*
// significance tests (check for overlapping confidence intervals)
capture log c
log using "T-tests.log", replace
foreach v of varlist `var' {
	mi estimate: mean `v', over(spmom)
}
log c
*/

preserve
mi extract 1, clear

// mean & sd by family structure
forvalues k = 0/1 {
	local i = 1
	foreach v of varlist `var' { 
		reg x `v' if spmom==`k', noconstant
		estadd summ
		mat max = e(max)
			if max[1,1] == 1 {
				mat d`k'`i' = nullmat(d`k'`i') \ e(mean), .
			}
			else if max[1,1] != 1 {
				mat d`k'`i' = nullmat(d`k'`i') \ e(mean), e(sd)
			}
		mat rowname d`k'`i' = `v'
		mat c`k' = nullmat(c`k')\ d`k'`i'
		local ++i
	}
}
restore

// min, max, & % imputed
preserve
mi extract 0, clear
local n = _N
foreach v of varlist `var' {
	reg x `v', noconstant
	estadd summ
	mat d = nullmat(d)\ e(min), e(max), ((`n'-e(N))/`n'*100)
}
fre spmom
restore
mat c = c0,c1,d

estout matrix(c, fmt(2 2 2 2 0 0 0)), nolz mlabel(none) collabels(Mean^a SD Mean^a SD Min Max "% Imputed") ///
	prehead("Table 1. Mean, Standard Deviation, Minimum, and Maximum of Study Variables, by Family Structure") ///
	varlabels(female "   Female" age "   Age (in years)" white "   ...White" black "   ...Black" asian "   ...Asian" natam "   ...Native American" hisp "   ...Hispanic" ///
			  nohs "   ...No HS degree" hsgrad "   ...HS degree" somecol "   ...Some college" colgrad "   ...College degree" inclog "   Income (logged)" ///
			  full "   Full-time employed parent" ///
			  h1id1f   "   ...We would hold hands" ///
			  h1id1l   "   ...We would talk about contraception or STD's" ///
			  h1id1m   "   ...We would kiss" ///
			  h1id1n   "   ...We would touch each other under our clothing or with no clothes on" ///
			  h1id1o   "   ...We would have sex" ///
			  h1id1p   "   ...My partner or I would get pregnant" ///
			  h1id1q   "   ...We would get married" ///
			  know     "Sexual knowledge quiz" ///
			  sex_phy  "   ...Physical pleasure" ///
			  sex_soc  "   ...Social promoters" ///
			  sex_resp "   ...Social inhibitors" ///
			  h1co1    "I have had sexual intercourse" ///
			  birthc1  "External barriers to birth control" ///
			  h1bc7    "Birth control is morally wrong" ///
			  h1bc8    "Friends will think I want sex if I use birth control" ///
			  birthc2  "Birth control self-efficacy^b" ///
			  h1co13   "   ...I have taken birth control pills regularly for at least one cycle" ///
			  h1co15   "   ...I currently take birth control pills" ///
			  h1co3    "My partner or I used birth control the first time I had sex^d" ///
			  h1co6    "My partner or I used birth control the last time I had sex^d" ///
			  h1rp2    "It wouldn't be bad to get pregnant" ///
			  prego    "Motivations for pregnancy" ///
			  sti      "A doctor has told me that I have an STD" ///
			  spmom	   "Single mother") ///	
	refcat(h1id1f "Ideal romantic relationship" sex_phy "Motivations for sexual intercourse" ///
	       h1co13 "Females only^c" female "Controls" white "Race-ethnicity" nohs "Parent educational attainment", label(" ")) ///
	postfoot("^a Proportions for dichotomous variables." ///
			 "^b Excludes 21 respondents who reported never wanting to use birth control." ///
			 "^c Includes only respondents who reported having sex (N = 1,443)." ///
			 "^d Includes only female respondents who reported having sex (N = 722)." ///
			 "Source: National Longitudinal Study of Adolescent to Adult Health.") ///
	sty(tab) delimit(;) // varwidth(72)




		////////////////////////////////////////////////////////////////////////////////////////////
		// program to merge models to report bivariate & multivariate results in their own column //
		////////////////////////////////////////////////////////////////////////////////////////////
		capt prog drop mergemodels								

		prog mergemodels, eclass
			version 8
			syntax namelist														// the program expects a list of names that identify models to merge
			tempname b V tmp													// set temporary names for matrices
			local i = 1															// local to use as a counter
			foreach name of local namelist {									// loop through each model listed in the list of names
				qui est restore `name'											// restore previously saved results to pull results
				local names "`names'`e(depvar)' "								// store the name of the dependent variable to use as matrix row/column names
				mat `b'`i' = e(b)												// put e(b) results in matrix to use
				mat `b' = nullmat(`b') , `b'`i'[1,1]							// put the coefficient in a vector to use in estout
				mat `tmp'`i' = e(V)												// put e(V) results in matrix to use
				capt confirm matrix `V'											// ask if matrix `V' exists
				if _rc {														// if matrix `V' doesn't exist...
					mat `V' = `tmp'`i'[1,1]										// 		...assign the variance to `V'
				}
				local j = `i'-1													// this local tells how many 0s need to go in `V' to add the new variance
				else {															// if matrix `V' exists...
					mat `V' = (`V',J(`j',1,0)) \ (J(1,`j',0) , `tmp'`i'[1,1] )	// 		...add 0s to the matrix with the new variance on the diagonal
				}
				local ++i
			}
		mat coln `b' = `names'													// these 3 name columns and rows as needed
		mat coln `V' = `names'
		mat rown `V' = `names'
		eret post `b' `V'														// this puts the matrices in e-class results to be used by estout
		eret local cmd "whatever"

		end
		////////////////////////////////////////////////////////////////////////////////////////////
		// program to merge models to report bivariate & multivariate results in their own column //
		////////////////////////////////////////////////////////////////////////////////////////////

capture log close
log using "Table 2 SEX.log", replace

local i = 1
foreach v of varlist h1id1f h1id1l h1id1m h1id1n h1id1o h1id1p h1id1q {
	display "`v'"
	eststo t2`i'a: mi estimate, post: logistic `v' spmom
		local a "`a't2`i'a "
	eststo t2`i'b: mi estimate, post: logistic `v' spmom female age ib1.race ib1.hpaed inclog full
		local b "`b't2`i'b "
	local ++i
}
foreach v of varlist know sex_phy sex_soc sex_resp {
	display "`v'"
	eststo t2`i'a: mi estimate, post: reg `v' spmom 
		local a "`a't2`i'a "
	eststo t2`i'b: mi estimate, post: reg `v' spmom female age ib1.race ib1.hpaed inclog full
		local b "`b't2`i'b "
	local ++i
}
	eststo t2`i'a: mi estimate, post: logistic h1co1 spmom
		local a "`a't2`i'a "
	eststo t2`i'b: mi estimate, post: logistic h1co1 spmom female age ib1.race ib1.hpaed inclog full
		local b "`b't2`i'b "
	local ++i
foreach v of varlist birthc1 h1bc7 h1bc8 birthc2 { 
	display "`v'"
	eststo t2`i'a: mi estimate, post: reg `v' spmom 
		local a "`a't2`i'a "
	eststo t2`i'b: mi estimate, post: reg `v' spmom female age ib1.race ib1.hpaed inclog full
		local b "`b't2`i'b "
	local ++i
}
preserve
	mi convert wide, clear
	keep if h1co1==1					// keeps only those who reported having sex
	foreach v of varlist h1co3 h1co6 {
		display "`v'"
		eststo t2`i'a: mi estimate, post: logistic `v' spmom 
			local a "`a't2`i'a "
		eststo t2`i'b: mi estimate, post: logistic `v' spmom female age ib1.race ib1.hpaed inclog full
			local b "`b't2`i'b "
		local ++i
	}
	foreach v of varlist h1co13 h1co15 {
		display "`v'"
		eststo t2`i'a: mi estimate, post: logistic `v' spmom if female==1
			local a "`a't2`i'a "
		eststo t2`i'b: mi estimate, post: logistic `v' spmom age ib1.race ib1.hpaed inclog full if female==1
			local b "`b't2`i'b "
		local ++i
	}
restore
foreach v of varlist h1rp2 prego { 
	display "`v'"
	eststo t2`i'a: mi estimate, post: reg `v' spmom 
		local a "`a't2`i'a "
	eststo t2`i'b: mi estimate, post: reg `v' spmom female age ib1.race ib1.hpaed inclog full
		local b "`b't2`i'b "
	local ++i
}
	eststo t2`i'a: mi estimate, post: logistic sti spmom
		local a "`a't2`i'a "
	eststo t2`i'b: mi estimate, post: logistic sti spmom female age ib1.race ib1.hpaed inclog full
		local b "`b't2`i'b "
	
local N = `e(N_mi)'
eststo bivariate:    mergemodels `a'
eststo multivariate: mergemodels `b'

reg birthc2 spmom 											if _mi_m==1
reg birthc2 spmom female age ib1.race ib1.hpaed inclog full if _mi_m==1

estout bivariate multivariate, cells("b(star fmt(3)) p" ci(par)) nolz ///
	mgroups("              Bivariate^a" "           Multivariate^a,b", pattern(1 1) span) mlabels(none) ///
	prehead("Table 2.""Single-parenthood and sex attitudes and behaviors") ///
	varlabels(h1id1f   "   ...We would hold hands" ///
			  h1id1l   "   ...We would talk about contraception or STD's" ///
			  h1id1m   "   ...We would kiss" ///
			  h1id1n   "   ...We would touch each other under our clothing or with no clothes on" ///
			  h1id1o   "   ...We would have sex" ///
			  h1id1p   "   ...My partner or I would get pregnant" ///
			  h1id1q   "   ...We would get married" ///
			  know     "Sexual knowledge quiz" ///
			  sex_phy  "   ...Physical pleasure" ///
			  sex_soc  "   ...Social promoters" ///
			  sex_resp "   ...Social inhibitors" ///
			  h1co1    "I have had sexual intercourse^c" ///
			  birthc1  "External barriers to birth control" ///
			  h1bc7    "Birth control is morally wrong" ///
			  h1bc8    "Friends will think I want sex if I use birth control" ///
			  birthc2  "Birth control self-efficacy^d" ///
			  h1co13   "   ...I have taken birth control pills regularly for at least one cycle" ///
			  h1co15   "   ...I currently take birth control pills" ///
			  h1co3    "My partner or I used birth control the first time I had sex^c,e" ///
			  h1co6    "My partner or I used birth control the last time I had sex^c,e" ///
			  h1rp2    "It wouldn't be bad to get pregnant" ///
			  prego    "Motivations for pregnancy" ///
			  sti      "A doctor has told me that I have an STD^c") ///
	refcat(h1id1f "Ideal romantic relationship^c" sex_phy "Motivations for sexual intercourse" ///
	       h1co13 "Females only^c,f", label(" ")) ///
	postfoot("Note:""Unstandardized coefficients from OLS regression unless otherwise indicated." ///
			 "N = `N' unless otherwise indicated." ///
			 "^a Single father is the reference category for structure (i.e., single moms = 1, single dads = 0)" ///
			 "^b Controls include gender, age, race, parent education, logged parental income, and parental employment." ///
			 "^c Log odds from logistic regression." ///
			 "^d Excludes 21 respondents who reported never wanting to use birth control." ///
			 "^e Includes only respondents who reported having sex (N = 1,443)." ///
			 "^f Includes only female respondents who reported having sex (N = 722)." ///
			 "Source: National Longitudinal Study of Adolescent to Adult Health.") ///
	sty(tab) delimit(;) // varwidth(72)
log close
		


/////////////////////////////
//                         //
//      interactions       //
//                         //
/////////////////////////////


		/////////////////////////////////////////////////////////////////////////////////////////
		// program to merge models to report main effects and interactions in their own column //
		/////////////////////////////////////////////////////////////////////////////////////////
		capt prog drop mergeint								

		prog mergeint, eclass
			version 8
			syntax namelist														// the program expects a list of names that identify models to merge
			tempname b V tmp													// set temporary names for matrices
			local i = 1															// local to use as a counter
			foreach name of local namelist {									// loop through each model listed in the list of names
				qui est restore `name'											// restore previously saved results to pull results
				local names "`names'spmom:`e(depvar)' female:`e(depvar)' int:`e(depvar)' "								// store the name of the dependent variable to use as matrix row/column names
				mat `b'`i' = e(b)												// put e(b) results in matrix to use
				mat `b' = nullmat(`b') , `b'`i'[1,1..3]							// put the coefficient in a vector to use in estout
				mat `tmp'`i' = e(V)												// put e(V) results in matrix to use
				capt confirm matrix `V'											// ask if matrix `V' exists
				if _rc {														// if matrix `V' doesn't exist...
					mat `V' = `tmp'`i'[1..3,1..3]								// 		...assign the variance to `V'
				}
				else {															// if matrix `V' exists...
					local j = 3*(`i'-1)											// this local tells how many 0s need to go in `V' to add the new variance
					mat `V' = (`V' , ///
							  J(`j',3,0)) \ (J(3,`j',0) , `tmp'`i'[1..3,1..3] )	// 		...add 0s to the matrix with the new variance on the diagonal
				}
				local ++i
			}
		mat coln `b' = `names'													// these 3 name columns and rows as needed
		mat coln `V' = `names'
		mat rown `V' = `names'
		eret post `b' `V'														// this puts the matrices in e-class results to be used by estout
		eret local cmd "whatever"

		end
		/////////////////////////////////////////////////////////////////////////////////////////
		// program to merge models to report main effects and interactions in their own column //
		/////////////////////////////////////////////////////////////////////////////////////////


capture log close
log using "Table 3 SEX.log", replace

local i = 1
foreach v of varlist h1id1f h1id1l h1id1m h1id1n h1id1o h1id1p h1id1q {
	display "`v'"
	eststo int`i': mi estimate, post: logistic `v' c.spmom##c.female age ib1.race ib1.hpaed inclog full
		local int "`int'int`i' "
	local ++i
}
foreach v of varlist know sex_phy sex_soc sex_resp {
	display "`v'"
	eststo int`i': mi estimate, post: reg `v' c.spmom##c.female age ib1.race ib1.hpaed inclog full
		local int "`int'int`i' "
	local ++i
}
	eststo int`i': mi estimate, post: logistic h1co1 c.spmom##c.female age ib1.race ib1.hpaed inclog full
		local int "`int'int`i' "
	local ++i
foreach v of varlist birthc1 h1bc7 h1bc8 birthc2 { 
	display "`v'"
	eststo int`i': mi estimate, post: reg `v' c.spmom##c.female age ib1.race ib1.hpaed inclog full
		local int "`int'int`i' "
	local ++i
}
preserve
	mi convert wide, clear
	keep if h1co1==1					// keeps only those who reported having sex
	foreach v of varlist h1co3 h1co6 {
		display "`v'"
		eststo int`i': mi estimate, post: logistic `v' c.spmom##c.female age ib1.race ib1.hpaed inclog full
		local int "`int'int`i' "
		local ++i
	}
	// interactions for h1co13 and h1co15 not included because they are only relevant for females anyway
restore
foreach v of varlist h1rp2 prego { 
	display "`v'"
	eststo int`i': mi estimate, post: reg `v' spmom female age ib1.race ib1.hpaed inclog full
		local int "`int'int`i' "
	local ++i
}
	eststo int`i': mi estimate, post: logistic sti c.spmom##c.female age ib1.race ib1.hpaed inclog full
		local int "`int'int`i' "
local N = e(N_mi)
eststo interact: mergeint `int'

estout interact, cells(b(fmt(3) star) ci(par fmt(3))) nolz unstack ///
	eqlabel("Single Mom" "Female" "Interaction") collabel(none) mlabel(none) ///
	prehead("Table 3. Single-parenthood and sex attitudes and behaviors: The interaction between parent and child sex") ///
	varlabels(h1id1f   "   ...We would hold hands" ///
			  h1id1l   "   ...We would talk about contraception or STD's" ///
			  h1id1m   "   ...We would kiss" ///
			  h1id1n   "   ...We would touch each other under our clothing or with no clothes on" ///
			  h1id1o   "   ...We would have sex" ///
			  h1id1p   "   ...My partner or I would get pregnant" ///
			  h1id1q   "   ...We would get married" ///
			  know     "Sexual knowledge quiz" ///
			  sex_phy  "   ...Physical pleasure" ///
			  sex_soc  "   ...Social promoters" ///
			  sex_resp "   ...Social inhibitors" ///
			  h1co1    "I have had sexual intercourse^a" ///
			  birthc1  "External barriers to birth control" ///
			  h1bc7    "Birth control is morally wrong" ///
			  h1bc8    "Friends will think I want sex if I use birth control" ///
			  birthc2  "Birth control self-efficacy^b" ///
			  h1co3    "My partner or I used birth control the first time I had sex^a,c" ///
			  h1co6    "My partner or I used birth control the last time I had sex^a,c" ///
			  h1co13   "   ...I have taken birth control pills regularly for at least one cycle" ///
			  h1co15   "   ...I currently take birth control pills" ///
			  h1rp2    "It wouldn't be bad to get pregnant" ///
			  prego    "Motivations for pregnancy" ///
			  sti      "A doctor has told me that I have an STD") ///
	refcat(h1id1f "Ideal romantic relationship^a" sex_phy "Motivations for sexual intercourse" ///
	       h1co13 "Females only^a,d", label(" ")) ///
	postfoot("Note:""Unstandardized coefficients from OLS regression unless otherwise indicated." ///
			 "N = `N' unless otherwise indicated." ///
			 "Controls include age, race, parent education, logged parental income, and parental employment." ///
			 "^a Log odds from logistic regression." ///
			 "^b Excludes 21 respondents who reported never wanting to use birth control." ///
			 "^c Includes only respondents who reported having sex (N = 1,443)." ///
			 "Source: National Longitudinal Study of Adolescent to Adult Health.") ///
	sty(tab) delimit(;) // varwidth(72)

log c


//
// supplementary
//

// motivations for sexual intercourse by item

// h1mo5 h1mo6 			are Physical
// h1mo1 h1mo7 h1mo8	are Social promoters
// h1mo2 h1mo3 h1mo4	are Social inhibitors
local a
local b
local i = 1

foreach v of varlist h1mo5 h1mo6 h1mo1 h1mo7 h1mo8 h1mo2-h1mo4 {
	display "`v'"
	eststo t2`i'a: mi estimate, post: reg `v' spmom 
		local a "`a't2`i'a "
	eststo t2`i'b: mi estimate, post: reg `v' spmom female age ib1.race ib1.hpaed inclog full
		local b "`b't2`i'b "
	local ++i
}

eststo bi:    mergemodels `a'
eststo multi: mergemodels `b'

estout bi multi, cells("b(star fmt(3)) p") nolz ///
	mgroups("              Bivariate^a" "           Multivariate^a,b", pattern(1 1) span) mlabels(none) ///
	prehead("Table 2 supplement.""Single-parenthood and motivations for having sex:""Unstandardized coefficients from OLS regression") ///
	varlabels(h1mo5   "...Would give me a great deal of pleasure" ///
			  h1mo6   "...It would relax me" ///
			  h1mo1   "...Friends would respect me more" ///
			  h1mo7   "...Would make me more attractive to opposite sex" ///
			  h1mo8   "...I would feel less lonely" ///
			  h1mo2   "...Partner would lose respect for me" ///
			  h1mo3   "...After sex I would feel guilty" ///
			  h1mo4   "...It would upset my mother") ///
	refcat(h1mo5 "Physical" h1mo1 "Social Promoters" h1mo2 "Social inhibitors", label(" ")) ///
	postfoot("Note:""Unstandardized coefficients from OLS regression unless otherwise indicated." ///
			 "N = 2,570" ///
			 "^a Single father is the reference category for structure (i.e., single moms = 1, single dads = 0)" ///
			 "^b Controls include gender, age, race, parent education, logged parental income, and parental employment" ///
			 "Source: National Longitudinal Study of Adolescent to Adult Health") ///
	/*sty(tab) delimit(;) */ varwidth(50)
