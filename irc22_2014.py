"""
Module for IRC 22:2014 bridge design clauses.
@author: Sweta Pal

"""

class IRC22_2014:

    # 601.4 Material Strength and Partial Safety Factor (Clause 601.4)

    @staticmethod
    def cl_601_4_material_safety_factors():
        """
        IRC:22-2014 Clause 601.4
        Table 1: Material Safety Factors (γm)

        Returns:
            dict containing partial safety factors for:
                - Steel yield stress
                - Steel ultimate stress
                - Reinforcement
                - Shear concrete
                - Bolts & Rivets
                - Welds
                - Concrete (accidental)
                (Separate factors for ultimate limit & serviceability)
        """

        # values as per Table 1 IRC-22:2014
        gamma_m = {
            "structural_steel_yield": {"ULS": 1.10, "SLS": 1.00},
            "structural_steel_ultimate": {"ULS": 1.25, "SLS": 1.00},
            "steel_reinforcement_yield": {"ULS": 1.15, "SLS": 1.00},
            "shear_concrete": {"ULS": 1.25, "SLS": 1.00},
            "bolts_rivets_shear_tension": {"ULS": 1.25, "SLS": 1.00},
            "welds_shop": {"ULS": 1.25, "SLS": 1.00},
            "welds_field": {"ULS": 1.50, "SLS": 1.00},
            "concrete_basic_seismic": {"ULS": 1.50, "SLS": 1.00},
            "concrete_accidental": {"ULS": 1.20, "SLS": 1.00}
        }

        return gamma_m

    # 602 Material  |  Annex III (Properties of Structural Steel)

    @staticmethod
    def cl_602_annexIII_structural_steel_properties():
        """
        IRC:22-2014
        Clause 602 | Annex-III

        III.1 Structural Steel:
        Properties to be assumed for all grades of steel for design.

        Returns:
            dict containing fundamental steel properties:
                - Young's Modulus
                - Shear Modulus
                - Poisson's Ratio
                - Coefficient of Thermal Expansion
        """

        properties = {
            "E_MPa": 2.0e5,                 # Young's Modulus (MPa)
            "E_GPa": 200,                   # Young's Modulus (GPa) - convenience
            "G_MPa": 0.77e5,                # Shear Modulus (MPa)
            "nu": 0.30,                     # Poisson's Ratio
            "thermal_expansion_per_C": 0.0000117  # per °C per unit length
        }

        return properties


    # 602 Material | Annex III (Concrete Properties)

    @staticmethod
    def cl_602_annexIII_concrete_properties(aggregate_type="quartzite"):
        """
        IRC:22-2014
        Clause 602 | Annex-III | III.2 Concrete

        Returns properties of concrete grades as per Table III-1:
            - fck : Cube compressive strength (MPa)
            - fcy : Cylinder compressive strength (MPa)
            - fctm: Mean tensile strength (MPa)
            - Ec  : Secant Modulus of Elasticity (GPa)

        Aggregate Modification (as per IRC note):
            quartzite/granite -> 1.0  (default)
            limestone          -> 0.9
            sandstone          -> 0.7
            basalt             -> 1.2
        """

        # Aggregate factors
        agg_factor_map = {
            "quartzite": 1.0,
            "granite": 1.0,
            "limestone": 0.9,
            "sandstone": 0.7,
            "basalt": 1.2
        }

        factor = agg_factor_map.get(aggregate_type.lower(), 1.0)

        # Table III-1 Data
        concrete = {
            "M15":  {"fck": 15, "fcy": 12, "fctm": 1.6, "Ec": 27},
            "M20":  {"fck": 20, "fcy": 16, "fctm": 1.9, "Ec": 29},
            "M25":  {"fck": 25, "fcy": 20, "fctm": 2.2, "Ec": 30},
            "M30":  {"fck": 30, "fcy": 24, "fctm": 2.8, "Ec": 31},
            "M35":  {"fck": 35, "fcy": 28, "fctm": 3.0, "Ec": 32},
            "M40":  {"fck": 40, "fcy": 32, "fctm": 3.3, "Ec": 33},
            "M45":  {"fck": 45, "fcy": 36, "fctm": 3.4, "Ec": 34},
            "M50":  {"fck": 50, "fcy": 40, "fctm": 3.5, "Ec": 35},
            "M55":  {"fck": 55, "fcy": 44, "fctm": 3.6, "Ec": 36},
            "M60":  {"fck": 60, "fcy": 48, "fctm": 3.7, "Ec": 37},
            "M65":  {"fck": 65, "fcy": 52, "fctm": 4.0, "Ec": 38},
            "M70":  {"fck": 70, "fcy": 56, "fctm": 4.4, "Ec": 38},
            "M75":  {"fck": 75, "fcy": 60, "fctm": 4.5, "Ec": 39},
            "M80":  {"fck": 80, "fcy": 64, "fctm": 4.7, "Ec": 40},
            "M85":  {"fck": 85, "fcy": 68, "fctm": 4.9, "Ec": 40},
            "M90":  {"fck": 90, "fcy": 72, "fctm": 5.0, "Ec": 41},
        }

        # Apply Aggregate Modulus Factor
        for grade, props in concrete.items():
            props["Ec"] = round(props["Ec"] * factor, 2)

        return concrete

    # 602 Material | Annex III | III.3 Reinforcement Steel

    @staticmethod
    def cl_602_annexIII_reinforcement_steel_properties():
        """
        IRC:22-2014
        Clause 602 | Annex-III | III.3 Reinforcement Steel

        Reinforcing steel used in concrete construction shall conform to:
            IS:1786 - 2008 (High Strength Deformed Bars & Wires)

        Table covered: Mechanical Properties of High Strength Deformed Bars

        Returns:
            dict of reinforcement grades with:
                fy  : characteristic yield strength (MPa)
                fu  : ultimate tensile strength (MPa)
                elongation_min_percent : minimum elongation %
                fu_by_fy_requirement  : minimum required tensile/yield ratio
        """

        reinforcement = {
            "Fe415": {
                "fy": 415,
                "fu": 485,
                "elongation_min_percent": 14.5,
                "fu_by_fy_requirement": 1.08
            },
            "Fe500": {
                "fy": 500,
                "fu": 565,
                "elongation_min_percent": 12.0,
                "fu_by_fy_requirement": 1.08
            },
            "Fe500D": {
                "fy": 500,
                "fu": 565,
                "elongation_min_percent": 16.0,
                "fu_by_fy_requirement": 1.08
            },
            "Fe550": {
                "fy": 550,
                "fu": 600,
                "elongation_min_percent": 10.0,
                "fu_by_fy_requirement": 1.08
            },
            "Fe550D": {
                "fy": 550,
                "fu": 600,
                "elongation_min_percent": 14.0,
                "fu_by_fy_requirement": 1.08
            },
            "Fe600": {
                "fy": 600,
                "fu": 660,
                "elongation_min_percent": 10.0,
                "fu_by_fy_requirement": 1.08
            },
        }

        return reinforcement

# table 2 
# to be written

    # 603.2 Effective Width of Concrete Slab
    # 603.2.1 Simply Supported Girder (Pinned at Both Ends)


    # 603.2.1  Effective Width of Simply Supported Girder

    @staticmethod
    def cl_603_2_1_effective_width_simply_supported(
        Lo,                  # Effective span = distance between zero moments (≈ L)
        beam_type="inner",   # "inner" or "outer"
        B1=None,             # spacing left of beam
        B2=None,             # spacing right of beam
        B0=None              # overhang / spacing outside outer girder
    ):
        """
        IRC:22-2014
        Clause 603.2.1 – Effective Width of Simply Supported Composite Girder

        Parameters:
            Lo  : Effective span (m or mm — maintain consistency)
            beam_type : "inner" or "outer"
            B1  : spacing to adjacent slab/girder on one side (for inner/outer)
            B2  : spacing to other side (for inner)
            B0  : spacing outside outer girder (for outer)

        Returns:
            dict {
                "beff": effective width,
                "clause": "IRC 22:2014 - 603.2.1",
                "equation": "3.2 / 3.3 / 3.4"
            }
        """

        beam_type = beam_type.lower()


        # INNER GIRDER  (Eq. 3.2 / 3.3)

        if beam_type == "inner":

            if B1 is None or B2 is None:
                raise ValueError("B1 and B2 must be provided for inner beams")

            # Eq 3.2
            beff = Lo / 4.0

            # Clause limit
            limit = (B1 + B2) / 2.0

            beff = min(beff, limit)

            return {
                "beam_type": "inner",
                "Lo": Lo,
                "B1": B1,
                "B2": B2,
                "beff": beff,
                "equation": "3.2 / 3.3",
                "clause": "IRC 22:2014 - 603.2.1"
            }


        # OUTER EDGE GIRDER (Eq. 3.4)

        elif beam_type == "outer":

            if B1 is None or B0 is None:
                raise ValueError("B1 and B0 must be provided for outer beams")

            X = B0

            # From clause: constraints
            # Lo/8 ≤ B1/2
            # B0 = X ≤ Lo/8
            if Lo / 8.0 > (B1 / 2.0):
                print("Warning: Condition Lo/8 ≤ B1/2 not satisfied as per IRC 22")
            if X > Lo / 8.0:
                print("Warning: X (=B0) exceeds Lo/8 limit")

            # Eq 3.4
            beff = Lo / 8.0 + X

            # upper bound:
            limit = (B1 / 2.0) + X
            beff = min(beff, limit)

            return {
                "beam_type": "outer",
                "Lo": Lo,
                "B1": B1,
                "B0": B0,
                "X": X,
                "beff": beff,
                "equation": "3.4",
                "clause": "IRC 22:2014 - 603.2.1"
            }

        else:
            raise ValueError("beam_type must be 'inner' or 'outer'")


    @staticmethod
    def cl_603_3_1_positive_moment_capacity(
            fck,
            fy,
            beff,
            ds,
            As,
            Af,
            bf,
            tf,
            tw,
            dc,
            combination_type="basic",
            is_compact=True,
            beff_compact_limit=None
        ):
            
            """
            IRC:22-2014
            Clause 603.3.1 + Annexure I (I.1)
            Positive Moment Resistance of Composite Beam
            Plastic / Compact Section with FULL Shear Interaction
            """

            # Annex I.2  —  Non-Compact Section Restriction
            # For non-compact sections, beff must be restricted
            # to the compact section limiting value from Table-2

            beff_used = beff   # default: full effective width

            if not is_compact:
                if beff_compact_limit is None:
                    raise ValueError(
                        "Non-compact section detected but beff_compact_limit not provided "
                        "(Table 2 compact limiting width must be supplied)."
                    )
                beff_used = min(beff, beff_compact_limit)


            # MATERIAL PARTIAL SAFETY FACTORS

            gamma_m = 1.10  # steel

            if combination_type in ["basic", "seismic"]:
                gamma_c = 1.50
            elif combination_type == "accidental":
                gamma_c = 1.20
            else:
                raise ValueError("Invalid load combination type")


            # CONCRETE PARAMETERS

            alpha_cc = 0.67

            if fck <= 60:
                eta = 1.0
                lam = 0.8
            else:
                eta = 1.0 - (fck - 60) / 250.0
                lam = 0.8 - (fck - 60) / 500.0


            # PARAMETER 'a'

            a = (alpha_cc * eta * fck / gamma_c) / (fy / gamma_m)


            # CASE SELECTION CONDITIONS

            left = beff_used * ds
            middle = a * As
            right = (bf * ds + 2 * a * Af)

            if left >= middle:
                case = 1     # NA in slab
            elif left < middle < right:
                case = 2     # NA in steel flange
            else:
                case = 3     # NA in web


            # xu COMPUTATION (Depth of NA)

            if case == 1:
                xu = (a * As) / beff_used

            elif case == 2:
                xu = ds + ((a * As - bf * ds) / (2 * bf * a))

            else:  # case 3
                xu = ds + tf + (a * (As - 2 * Af) - bf * ds) / (2 * a * tw)


            # ULTIMATE MOMENT CAPACITY Mp

            fy_eff = fy / gamma_m

            if case == 1:
                Mp = As * fy_eff * (dc + 0.5 * ds - lam * xu / 2)

            elif case == 2:
                Mp = fy_eff * (
                    As * (dc + 0.5 * ds * (1 - lam))
                    - bf * (xu - ds) * tf
                    + (1 - lam) * ds * tf
                )

            else:  # case 3
                Mp = fy_eff * (
                    As * (dc + 0.5 * ds * (1 - lam))
                    - 2 * Af * (0.5 * tf * (1 - lam) * ds)
                    - tw * (xu - ds - tf) * tf
                )

            return {
                "gamma_m": gamma_m,
                "gamma_c": gamma_c,
                "alpha_cc": round(alpha_cc, 4),
                "eta": round(eta, 4),
                "lambda": round(lam, 4),
                "a": round(a, 4),
                "case": case,
                "xu_m": round(xu, 6),
                "Mp_kNm": round(Mp / 1e6, 3),  # assume Nmm → convert to kNm if input mm/N
                "clause": "IRC 22:2014 - 603.3.1 + Annex I (I.1 & I.2)"
            }

    @staticmethod
    def cl_603_3_3_1_buckling_resistance_moment(
        Mp,             # plastic moment from 603.3.1 (I.1 / I.2)
        fy,             # MPa
        Zp,             # plastic modulus (mm3)
        Z=None,         # elastic modulus (mm3) (needed for semi-compact)
        section_type="plastic",   # "plastic" / "compact" / "semi-compact"
        fabrication="rolled",      # "rolled" / "welded"
        E=2.0e5,        # MPa
        G=0.77e5,       # MPa
        Iy=None,        # mm4
        It=None,        # mm4
        Iw=None,        # mm6
        LLT=None        # mm (effective buckling length)
    ):
        """
        IRC:22-2014
        Clause 603.3.3.1
        Annexure I: I.5 Buckling Resistance Moment (Construction Stage)

        Returns:
            χLT and buckling reduced plastic moment
        """

        if None in [Iy, It, Iw, LLT]:
            raise ValueError("Iy, It, Iw and LLT must be provided")

        # βz

        if section_type.lower() in ["plastic", "compact"]:
            beta_z = 1.0
        else:   # semi compact
            if Z is None:
                raise ValueError("Z required for semi-compact sections")
            beta_z = Z / Zp

        # αLT

        if fabrication.lower() == "rolled":
            alpha_LT = 0.21
        else:
            alpha_LT = 0.49


        # Critical buckling moment Mcr

        import math

        term1 = (math.pi**2 * E * Iy) / (LLT**2)
        term2 = (G * It) + ((math.pi**2 * E * Iw) / (LLT**2))

        Mcr = math.sqrt(term1 * term2)   # Nmm

        # Slenderness ratio λLT
        lambda_LT = math.sqrt(beta_z * (Zp * fy) / Mcr)

        # If λLT <= 0.4 → no LTB effect

        if lambda_LT <= 0.4:
            chi_LT = 1.0
        else:
            phi_LT = 0.5 * (1 + alpha_LT * (lambda_LT - 0.2) + lambda_LT**2)
            chi_LT = 1 / (phi_LT + math.sqrt(phi_LT**2 - lambda_LT**2))
            chi_LT = min(chi_LT, 1.0)


        # Buckling reduced plastic moment
        
        Mpl_buck = chi_LT * Mp

        return {
            "beta_z": round(beta_z, 4),
            "alpha_LT": alpha_LT,
            "lambda_LT": round(lambda_LT, 4),
            "chi_LT": round(chi_LT, 4),
            "Mcr_kNm": round(Mcr / 1e6, 3),
            "Mpl_buckling_kNm": round(Mpl_buck / 1e6, 3),
            "clause": "IRC 22:2014 - 603.3.3.1 Annex I (I.5)"
        }


    # 603.3.3.2  VERTICAL SHEAR :  (1) PLASTIC SHEAR RESISTANCE

    @staticmethod
    def cl_603_3_3_2_plastic_shear_resistance(
        section_type,
        fyw,
        h=None,
        d=None,
        tw=None,
        A=None,
        b=None
    ):
        """
        IRC:22-2014
        Clause 603.3.3.2  (1) Plastic Shear Resistance
        Equation 3.5

        Vn = Vp = (Av * fyw) / √3
        Vd = Vn / γm0
        γm0 = 1.10

        Parameters:
            section_type : "I_major", "I_minor", "RHS", "CHS"
            fyw : yield stress of web (MPa)
            h   : total section depth (mm)           [for I hot rolled major]
            d   : clear web depth (mm)               [for welded I or minor]
            tw  : web thickness (mm)
            A   : total cross section area (mm2)     [for RHS / CHS]
            b   : flange width / tube breadth (mm)

        Returns:
            dict with:
                Av          : shear area (mm2)
                Vn_kN       : nominal shear resistance (kN)
                Vd_kN       : design shear resistance (kN)
                gamma_m0    : partial safety factor
                clause      : reference
        """

        import math

        gamma_m0 = 1.10

        section_type = section_type.lower()
        # Compute Shear Area Av

        if section_type == "i_major":
            # I / Channel – Major Axis
            # Av = h * tw  (hot rolled)
            # or d * tw    (welded). We assume user gives correct value
            if h is None or tw is None:
                raise ValueError("h and tw required for I_major section")
            Av = h * tw

        elif section_type == "i_minor":
            # Minor axis bending
            # Av = 2 * flange_width * flange_thickness
            if b is None or tw is None:
                raise ValueError("b and tw required for I_minor section")
            Av = 2 * b * tw

        elif section_type == "rhs":
            # Rectangular Hollow Section
            # Av = A * d / (b + d)
            if A is None or b is None or d is None:
                raise ValueError("A, b, d required for RHS")
            Av = A * d / (b + d)

        elif section_type == "chs":
            # Circular Hollow Section
            # Av = 2A / π
            if A is None:
                raise ValueError("A required for CHS")
            Av = (2 * A) / math.pi

        else:
            raise ValueError("Invalid section type")


        # Nominal Shear Strength

        Vn = (Av * fyw) / math.sqrt(3)

        # Design Shear Strength
        Vd = Vn / gamma_m0

        return {
            "Av_mm2": Av,
            "Vn_kN": round(Vn / 1000, 3),  # convert N to kN
            "Vd_kN": round(Vd / 1000, 3),
            "gamma_m0": gamma_m0,
            "clause": "IRC 22:2014 - 603.3.3.2 (1) Plastic Shear Resistance"
        }


    # 603.3.3.2 (2)(a)
    # Shear Buckling Resistance – Simple Post-Critical Method

    @staticmethod
    def cl_603_3_3_2_shear_buckling_post_critical(
        Av,          # shear area (mm2)
        fyw,         # web yield strength MPa
        d,           # web depth (mm)
        tw,          # web thickness (mm)
        c=None,      # spacing of transverse stiffeners (mm)
        E=2.0e5,     # MPa
        mu=0.3,      # Poisson ratio
        stiffeners_at_support_only=True
    ):
        """
        IRC:22-2014
        Clause 603.3.3.2 (2)(a)
        Simple Post-Critical Shear Buckling Method

        Returns:
            Nominal shear strength Vn (kN)
        """

        import math

        # Step-1: kv Based on IRC Rules

        if stiffeners_at_support_only:
            kv = 5.35
        else:
            if c is None:
                raise ValueError("c (stiffener spacing) required when intermediate stiffeners exist")

            ratio = c / d

            if ratio < 1.0:
                kv = 4.0 + 5.35 / (ratio**2)
            else:
                kv = 5.35 + 4.0 / (ratio**2)


        # Step-2: Elastic Critical Shear Stress τcr,e
        tau_cr_e = (
            (kv * (math.pi**2) * E)
            / (12 * (1 - mu**2))
        ) * (tw / d) ** 2


        # Step-3: λw (Web Slenderness Ratio)

        lambda_w = math.sqrt(fyw / (math.sqrt(3) * tau_cr_e))

        # Step-4: τb Based on λw Regions
        tau_y = fyw / math.sqrt(3)

        if lambda_w <= 0.8:
            tau_b = tau_y

        elif lambda_w < 1.2:
            tau_b = (1 - 0.8 * (lambda_w - 0.8)) * tau_y

        else:
            tau_b = tau_y / (lambda_w**2)

        # Step-5: Final Nominal Shear Strength
        Vn = Av * tau_b   # N

        return {
            "kv": round(kv, 4),
            "tau_cr_e_MPa": round(tau_cr_e, 4),
            "lambda_w": round(lambda_w, 4),
            "tau_b_MPa": round(tau_b, 4),
            "Vn_kN": round(Vn / 1000, 3),
            "clause": "IRC 22:2014 - 603.3.3.2 (2)(a)"
        }

    @staticmethod
    def cl_603_3_3_2_tension_field_method(
            Av,          # shear area of web (mm2)
            tw,          # web thickness (mm)
            d,           # clear web depth between flanges (mm)
            c,           # spacing of transverse stiffeners (mm)
            fyw,         # yield stress of web (MPa)
            fyf,         # yield stress of flange (MPa)
            bf,          # flange width (mm)
            tf,          # flange thickness (mm)
            tau_b,       # buckling shear stress from 603.3.3.2 (2)(a)
            gamma_m0=1.10
    ):
        """
        IRC 22:2014
        Clause 603.3.3.2 (2)(b)
        Tension Field Method (IS 800 8.4.2.2(b))

        Conditions:
            - Transverse stiffeners at supports
            - Intermediate stiffeners present
            - Panels adjacent provide anchorage
            - c/d ≥ 1.0
        """

        import math

        # CHECK CONDITION 
        if c / d < 1.0:
            raise ValueError("Tension Field Method not permitted because c/d < 1.0")

        # STEP 1 : φ (inclination)
        phi = math.atan(d / c)   # radians

        #  STEP 2 : ψ parameter 
        psi = 1.5 * tau_b * math.sin(2 * phi)

        #  STEP 3 : Tension field yield strength fv 
        fv = math.sqrt(fyw**2 - 3 * (tau_b**2) + psi**2) - psi

        #  STEP 4 : Anchorage length s 
        # Mfr from Eq 3.12
        Mfr = 0.25 * bf * (tf**2) * fyf * (1 - (0 / (bf * tf * fyf / gamma_m0))**2)  # assuming no axial force Nf = 0

        s = (2 / math.sin(phi)) * math.sqrt(Mfr / (fyw * tw))
        s = min(s, c)   # must not exceed c

        #  STEP 5 : Width of tension field 
        wtf = d * math.cos(phi) + (c - s) * math.sin(phi)

        #  STEP 6 : Nominal shear strength Vtf 
        Vtf = Av * tau_b + 0.9 * wtf * tw * fv * math.sin(phi)

        #  STEP 7 : Plastic shear cap
        Vp = Av * fyw / math.sqrt(3)

        # governing
        Vn = min(Vtf, Vp)

        # DESIGN STRENGTH
        Vd = Vn / gamma_m0

        return {
            "phi_rad": round(phi, 6),
            "psi_MPa": round(psi, 4),
            "fv_MPa": round(fv, 4),
            "wtf_mm": round(wtf, 3),
            "anchorage_s_mm": round(s, 3),
            "Vtf_kN": round(Vtf / 1000, 3),
            "Vp_kN": round(Vp / 1000, 3),
            "Vn_kN": round(Vn / 1000, 3),
            "Vd_kN": round(Vd / 1000, 3),
            "clause": "IRC 22:2014 - 603.3.3.2 (2)(b) Tension Field Method"
        }
    

    @staticmethod
    def cl_603_3_3_3_reduction_bending_high_shear(
        V,          # applied factored shear force (kN or N consistent)
        Vd,         # design shear strength (same unit as V)
        Md,         # plastic design moment of full section (Nmm or consistent)
        section_type="plastic",   # "plastic" / "compact" / "semi-compact"
        Mfd=None,   # plastic design moment excluding shear area (needed for plastic/compact)
        Ze=None,    # elastic section modulus (for semi-compact)
        fy=None,    # MPa
        gamma_m0=1.10
    ):
        """
        IRC 22:2014
        Clause 603.3.3.3 - Reduction in Bending Resistance Under High Shear

        Handles:
            - Plastic / Compact Section   → Eq. 3.13
            - Semi-Compact Section        → Eq. 3.14
        """

        section_type = section_type.lower()

        # -----------------------------
        # If V <= 0.6 Vd → No reduction
        # -----------------------------
        if V <= 0.6 * Vd:
            return {
                "reduction_required": False,
                "Mdv": Md,
                "Mdv_kNm": round(Md / 1e6, 3),
                "clause": "IRC 22:2014 - 603.3.3.3 (No reduction, V <= 0.6Vd)"
            }

        # ----------------------------------------
        # Case 1: Plastic or Compact Sections
        # ----------------------------------------
        if section_type in ["plastic", "compact"]:

            if Mfd is None:
                raise ValueError("Mfd must be provided for plastic/compact sections")

            beta = (2 * (V / Vd) - 1) ** 2

            Mdv = Md - beta * (Md - Mfd)

            # Upper limit 1.2 Ze fy / γm0
            if Ze is not None and fy is not None:
                upper_limit = 1.2 * Ze * fy / gamma_m0
                Mdv = min(Mdv, upper_limit)
            else:
                upper_limit = None

            return {
                "reduction_required": True,
                "beta": round(beta, 4),
                "Mdv": Mdv,
                "Mdv_kNm": round(Mdv / 1e6, 3),
                "upper_limit_kNm": None if upper_limit is None else round(upper_limit / 1e6, 3),
                "clause": "IRC 22:2014 - 603.3.3.3 (Eq. 3.13 Plastic/Compact)"
            }

        # ----------------------------------------
        # Case 2: Semi-Compact Section
        # ----------------------------------------
        elif section_type == "semi-compact":

            if Ze is None or fy is None:
                raise ValueError("Ze and fy must be provided for semi-compact section")

            Mdv = Ze * fy / gamma_m0

            return {
                "reduction_required": True,
                "Mdv": Mdv,
                "Mdv_kNm": round(Mdv / 1e6, 3),
                "clause": "IRC 22:2014 - 603.3.3.3 (Eq. 3.14 Semi-Compact)"
            }

        else:
            raise ValueError("section_type must be 'plastic', 'compact', or 'semi-compact'")

    
    # 604.3 Stresses and Deflection — Modular Ratio

    @staticmethod
    def cl_604_3_modular_ratio(
            Es=2.0e5,          # MPa (Steel Modulus)
            Ecm=None,          # MPa (Concrete modulus at 28 days)
            Eci=None,          # MPa (Concrete modulus before 28 days)
            Kc=0.5             # creep factor
    ):
        """
        IRC:22-2014
        Clause 604.3 - Stresses and Deflection
        Modular Ratio (m)

        Returns modular ratio:
            - Before 28 days (construction stage)
            - Short-term after setting
            - Long-term after setting
        """

        if Ecm is None:
            raise ValueError("Ecm (Concrete modulus at 28 days) must be provided")

        # ---------- BEFORE 28 DAYS ----------
        m_before = None
        if Eci:
            m_before = Es / Eci

        # ---------- AFTER 28 DAYS ----------
        # Short-term
        m_short = Es / Ecm
        m_short = max(m_short, 7.5)

        # Long-term
        m_long = Es / (Kc * Ecm)
        m_long = max(m_long, 15.0)

        return {
            "m_before_28days": m_before,
            "m_short_term": round(m_short, 3),
            "m_long_term": round(m_long, 3),
            "Es_MPa": Es,
            "Ecm_MPa": Ecm,
            "Eci_MPa": Eci,
            "Kc": Kc,
            "clause": "IRC 22:2014 - 604.3"
        }


    # 604.3.1  Limiting Stresses for Serviceability

    @staticmethod
    def cl_604_3_1_limiting_stresses(
        fck_concrete,
        fy_steel,
        fbc=None,   # bending compressive stress in steel (MPa)
        fbt=None,   # bending tensile stress in steel (MPa)
        fp=0,       # bearing stress in steel (MPa)
        tau_b=0     # shear stress in steel (MPa)
    ):
        """
        IRC:22-2014
        Clause 604.3.1 - Limiting Stresses for Serviceability

        Covers:
        1. Concrete compressive stress limit  -> IRC 112-2011 Cl. 12.2.1
        2. Reinforcement tensile stress limit -> IRC 112-2011 Cl. 12.2.2
        3. Structural steel equivalent stress -> must be <= 0.9 fy

        Parameters:
            fck_concrete (float): concrete characteristic cube strength (MPa)
            fy_steel (float): steel yield stress (MPa)

            fbc (float): actual compressive bending stress in steel (MPa)
            fbt (float): actual tensile bending stress in steel (MPa)
            fp (float):   bearing stress in steel section (MPa)
            tau_b (float): shear stress in steel section (MPa)

        Returns:
            dict containing allowable limits + equivalent stress check
        """

        import math


        # Concrete Limit (IRC 112-2011 Cl 12.2.1)
        k1 = 0.48
        f_conc_allow = k1 * fck_concrete


        # Reinforcement Steel Limit (IRC 112-2011 Cl 12.2.2)

        k3 = 0.80
        f_reinf_allow = k3 * fy_steel


        #  Structural Steel Equivalent Stress Limit
        steel_limit = 0.9 * fy_steel

        fe_comp = None
        fe_tens = None

        if fbc is not None and fbt is not None:

            # Equivalent compressive stress (En 4.1)
            fe_comp = math.sqrt(
                fbc**2 + fp**2 + fbc*fp + 3*(tau_b**2)
            )

            # Equivalent tensile stress (En 4.2)
            fe_tens = math.sqrt(
                fbt**2 + fp**2 + fbt*fp + 3*(tau_b**2)
            )

            steel_safe = (fe_comp <= steel_limit) and (fe_tens <= steel_limit)

        else:
            steel_safe = "Not evaluated (fbc/fbt not provided)"


        return {
            "concrete_allowable_stress_MPa": round(f_conc_allow, 3),
            "reinforcement_allowable_stress_MPa": round(f_reinf_allow, 3),
            "steel_equivalent_limit_0.9fy_MPa": round(steel_limit, 3),

            "steel_equivalent_compressive_fe_MPa": None if fe_comp is None else round(fe_comp, 3),
            "steel_equivalent_tensile_fe_MPa": None if fe_tens is None else round(fe_tens, 3),

            "steel_serviceability_safe": steel_safe,

            "clause": "IRC 22:2014 - 604.3.1 + IRC 112-2011 12.2.1 / 12.2.2"
        }


    # 604.3.2  Limit of Deflection and Camber
    @staticmethod
    def cl_604_3_2_deflection_limits(
        span_m,
        defl_live_impact_mm=None,
        defl_total_mm=None,
        cantilever_m=0,
        cant_defl_total_mm=None,
        cant_defl_live_impact_mm=None
    ):

        """
        IRC 22:2014
        Clause 604.3.2 - Limit of Deflection and Camber

        Checks:
        Main Span
            Live load + Impact  <= L / 800
            Total deflection    <= L / 600

        Cantilever Tip (if exists)
            Total deflection            <= Lc / 300
            Live load + Impact defl     <= Lc / 400

        Parameters:
            span_m (float): span of girder in meters
            defl_live_impact_mm (float): deflection under live + impact (mm)
            defl_total_mm (float): total deflection DL + SDL + LL + IMP (mm)

            cantilever_m (float): cantilever length (m), 0 if none
            cant_defl_total_mm (float): total tip deflection (mm)
            cant_defl_live_impact_mm (float): tip deflection LL + impact (mm)

        Returns:
            dict describing limits and pass / fail status
        """

        span_mm = span_m * 1000

        # MAIN GIRDER LIMITS
        allow_live_impact = span_mm / 800.0
        allow_total = span_mm / 600.0

        live_ok = None
        total_ok = None

        if defl_live_impact_mm is not None:
            live_ok = defl_live_impact_mm <= allow_live_impact

        if defl_total_mm is not None:
            total_ok = defl_total_mm <= allow_total


        # CANTILEVER CHECKS
        cantilever_results = None

        if cantilever_m > 0:
            Lc_mm = cantilever_m * 1000

            allow_cant_total = Lc_mm / 300.0
            allow_cant_live = Lc_mm / 400.0

            cant_total_ok = None
            cant_live_ok = None

            if cant_defl_total_mm is not None:
                cant_total_ok = cant_defl_total_mm <= allow_cant_total

            if cant_defl_live_impact_mm is not None:
                cant_live_ok = cant_defl_live_impact_mm <= allow_cant_live

            cantilever_results = {
                "cantilever_length_mm": Lc_mm,
                "allow_total_mm": round(allow_cant_total, 3),
                "allow_live_impact_mm": round(allow_cant_live, 3),
                "total_ok": cant_total_ok,
                "live_impact_ok": cant_live_ok
            }

        return {
            "span_mm": span_mm,

            "main_girder_limits": {
                "allow_live_impact_mm": round(allow_live_impact, 3),
                "allow_total_mm": round(allow_total, 3),
                "live_check_ok": live_ok,
                "total_check_ok": total_ok
            },

            "cantilever_check": cantilever_results,

            "clause": "IRC 22:2014 - 604.3.2"
        }

    @staticmethod
    def cl_604_4_crack_control_As_min(
        fctm,              # MPa (mean tensile strength of concrete)
        beff,              # mm effective flange width in tension
        t_slab,            # mm slab thickness in tension zone
        fy,                # MPa reinforcement yield
        kc=0.5,            # coefficient ≥ 0.5 (default recommended)
        k=1.0,             # restraint coefficient (1.0 typical)
        sigma_s=None,      # MPa steel stress to use, default = fy
        As_provided=None   # mm2 optional provided reinforcement to check adequacy
    ):
        """
        IRC 22:2014
        Clause 604.4 - Crack Control
        Ref: IRC:112-2011 Clause 12.3.3

        Computes minimum reinforcement required in tension flange region
        to control cracking.

        Returns:
            dict results including As_min and pass/fail (if As_provided given)
        """

        # Basic validity checks 
        if kc < 0.5:
            raise ValueError("kc must be ≥ 0.5 as per IRC:112-2011 guidance")

        if sigma_s is None:
            sigma_s = fy    # assume steel stress = fy unless user supplies lower value

        # Tensile Concrete Area 
        Act = beff * t_slab   # mm2

        # Minimum Reinforcement
        # As_min = kc * k * fct_eff * Act / sigma_s
        As_min = kc * k * fctm * Act / sigma_s

        result = {
            "Act_mm2": round(Act, 2),
            "kc": kc,
            "k": k,
            "fctm_MPa": fctm,
            "sigma_s_MPa": sigma_s,
            "As_min_mm2": round(As_min, 2),
            "clause": "IRC 22:2014 - 604.4 | IRC 112-2011 Cl.12.3.3"
        }

        # Optional adequacy check
        if As_provided is not None:
            result["As_provided_mm2"] = As_provided
            result["is_ok"] = As_provided >= As_min

        return result

    @staticmethod
    def cl_605_2_fatigue_design(
            tp,                 # thicker plate thickness (mm)
            f,                  # actual fatigue stress range (MPa)
            Nsc=None,           # actual no. of expected stress cycles
            gamma_mf=1.35,      # partial safety factor for fatigue strength
            gamma_m=1.0         # partial safety factor for loads
        ):
        """
        IRC:22-2015
        Clause 605.2 - Fatigue Design

        Handles:
        1) Capacity Reduction Factor for thickness > 25mm
        μr = (25/tp)^0.25 <= 1.0

        2) Low Fatigue Check:
        Fatigue assessment NOT REQUIRED if:
            f < 27 / γ_mf
        OR
            Nsc < 5e6 * ((27/γ_mf) / (γ_m * f))^3

        Parameters:
            tp          plate thickness (mm)
            f           fatigue stress range (MPa)
            Nsc         expected stress cycles
            gamma_mf    fatigue strength safety factor (default 1.35 typical worst)
            gamma_m     load factor (default 1.0)

        Returns dict with:
            mu_r
            stress_limit
            stress_ok
            cycle_limit
            cycle_ok
            fatigue_required (True/False)
        """

        import math

        # Capacity Reduction Factor μr

        if tp <= 25:
            mu_r = 1.0
        else:
            mu_r = (25.0 / tp) ** 0.25
            mu_r = min(mu_r, 1.0)

        # Low Fatigue - Stress Range Criterion
        stress_limit = 27.0 / gamma_mf
        stress_ok = f < stress_limit

        # Low Fatigue - Number of Cycles Criterion

        cycle_limit = None
        cycle_ok = None

        if Nsc is not None and f > 0:
            cycle_limit = 5e6 * ((27.0 / gamma_mf) / (gamma_m * f))**3
            cycle_ok = Nsc < cycle_limit

        # Final Decision
        # Fatigue NOT required if either condition passes
        if stress_ok or (cycle_ok is True):
            fatigue_required = False
        else:
            fatigue_required = True

        return {
            "mu_r": round(mu_r, 4),

        "stress_limit_27_by_gamma_mf": round(stress_limit, 3),
        "f_stress_range": f,
        "stress_condition_ok": stress_ok,

        "cycle_limit_Nsc_allowed": None if cycle_limit is None else round(cycle_limit),
        "Nsc_input": Nsc,
        "cycle_condition_ok": cycle_ok,

        "fatigue_required": fatigue_required,
        "clause": "IRC 22:2015 - 605.2 Fatigue Design"
    }

    @staticmethod
    def cl_605_3_fatigue_strength(
            Nsc,            # number of stress cycles
            ffn=None,       # normal fatigue strength at 5e6 cycles (MPa)
            tfn=None        # shear fatigue strength at 5e6 cycles (MPa)
        ):
        """
        IRC 22:2015
        Clause 605.3 - Fatigue Strength

        Computes:
            - Design normal fatigue stress range f_f
            - Design shear fatigue stress range tau_f

        Equations:
            For Normal Stress:
                If Nsc <= 5e6      →  f_f = ffn
                If 5e6 < Nsc ≤ 1e8 →  f_f = ffn * (5e6/Nsc)^(1/3)

            For Shear Stress:
                tau_f = tfn * (5e6/Nsc)^(1/5)

        Returns dictionary
        """

        import math

        if Nsc <= 0:
            raise ValueError("Nsc must be positive")

        # ---------- Normal Fatigue Strength ----------
        f_f = None
        if ffn is not None:
            if Nsc <= 5e6:
                f_f = ffn
            elif Nsc <= 1e8:
                f_f = ffn * (5e6 / Nsc) ** (1.0 / 3.0)
            else:
                # IRC usually implies curve decay beyond 1e8 as continuing trend
                f_f = ffn * (5e6 / Nsc) ** (1.0 / 3.0)

        # ---------- Shear Fatigue Strength ----------
        tau_f = None
        if tfn is not None:
            tau_f = tfn * (5e6 / Nsc) ** (1.0 / 5.0)

        return {
            "Nsc": Nsc,
            "f_f_normal_MPa": None if f_f is None else round(f_f, 3),
            "tau_f_shear_MPa": None if tau_f is None else round(tau_f, 3),
            "clause": "IRC 22:2015 - 605.3 Fatigue Strength"
        }


    @staticmethod
    def cl_605_4_fatigue_assessment(
        ff,          # Normal fatigue strength range for NSC (from 605.3)
        tf,          # Shear fatigue strength range for NSC (from 605.3)
        mu_r=1.0,    # Capacity reduction factor (605.2), default 1.0
        gamma_mft=1.35,  # Partial safety factor for fatigue strength (Table 3)
        fy=None,     # Yield stress of steel (MPa)
        f_actual=None,   # Actual normal stress range in service (MPa)
        tau_actual=None  # Actual shear stress range in service (MPa)
    ):
        """
        IRC:22-2015
        Clause 605.4 – Fatigue Assessment

        Returns design fatigue strength ranges and compliance checks.
        """

        if fy is None:
            raise ValueError("fy (yield stress) must be provided")

        # Design Fatigue Strengths
        f_fd = mu_r * ff / gamma_mft    # Eq 5.4
        tau_fd = mu_r * tf / gamma_mft  # Eq 5.5

        # Elastic Upper Limits
        f_limit_elastic = 1.5 * fy
        tau_limit_elastic = 1.5 * fy / (3 ** 0.5)

        results = {
            "f_fd_MPa": round(f_fd, 3),
            "tau_fd_MPa": round(tau_fd, 3),
            "elastic_limit_normal_MPa": round(f_limit_elastic, 3),
            "elastic_limit_shear_MPa": round(tau_limit_elastic, 3),
            "clause": "IRC 22:2015 - 605.4 Fatigue Assessment"
        }

        # If actual stress values are supplied, perform checks
        if f_actual is not None:
            results["normal_within_f_fd"] = f_actual <= f_fd
            results["normal_within_elastic_limit"] = f_actual <= f_limit_elastic

        if tau_actual is not None:
            results["shear_within_tau_fd"] = tau_actual <= tau_fd
            results["shear_within_elastic_limit"] = tau_actual <= tau_limit_elastic

        return results

    @staticmethod
    def cl_606_3_1_stud_connector_strength(
    d_mm,
    hs_mm,
    fu_MPa,
    fck_MPa,
    Ecm_MPa,
    gamma_v=1.25
    ):
        """
        IRC:22-2015
        Clause 606.3.1  – Ultimate strength of shear connectors (Stud type)

        Returns design resistance Qu (N)
        """

        import math

        # --- basic checks ---
        if d_mm <= 0 or hs_mm <= 0:
            raise ValueError("Stud diameter and height must be positive")

        if fu_MPa > 500:
            # as per clause <= 500 MPa recommended
            print("Warning: fu exceeds 500 MPa, IRC recommends fu ≤ 500 MPa")

        # cylinder strength of concrete 
        fck_cyl = 0.8 * fck_MPa

        # alpha factor 
        slenderness = hs_mm / d_mm

        if slenderness >= 4:
            alpha = 1.0
        elif 3 < slenderness < 4:
            alpha = 0.2 * (slenderness + 1)
        else:
            # below code range – normally not allowed
            print("Warning: hs/d < 3, outside IRC recommended range")
            alpha = 0.2 * (slenderness + 1)

        # steel governed capacity 
        Qu_steel = (0.8 * fu_MPa * math.pi * (d_mm**2) / 4.0) / gamma_v

        # concrete governed capacity 
        Qu_conc = (0.29 * alpha * (d_mm**2) * math.sqrt(fck_cyl * Ecm_MPa)) / gamma_v

        # governing value 
        Qu = min(Qu_steel, Qu_conc)

        return {
            "alpha": round(alpha, 4),
            "slenderness_h_by_d": round(slenderness, 3),
            "Qu_steel_N": round(Qu_steel, 2),
            "Qu_concrete_N": round(Qu_conc, 2),
            "Qu_governing_N": round(Qu, 2),
            "Qu_governing_kN": round(Qu / 1000.0, 3),
            "limit_state": "min(steel, concrete)",
            "clause": "IRC 22:2015 - 606.3.1 Stud Connectors"
        }

    @staticmethod
    def cl_606_3_2_shear_connector_fatigue_strength(
            Nsc,
            tau_fn=67.0
        ):
        """
        IRC:22-2015
        Clause 606.3.2 - Fatigue strength of shear connectors (Stud)
        
        τf = τfn * (5e6 / Nsc) ^ (1/5)

        Parameters
        ----------
        Nsc : float
            Number of stress cycles
        tau_fn : float
            Nominal fatigue shear strength at 5×10^6 cycles (MPa / N/mm2)
            Default = 67 MPa for stud connectors (Table 8)

        Returns
        -------
        dict with τf MPa
        """

        import math

        if Nsc <= 0:
            raise ValueError("Nsc must be positive")

        tau_f = tau_fn * ((5e6 / Nsc) ** 0.2)

        return {
            "tau_fn_MPa": tau_fn,
            "Nsc": Nsc,
            "tau_f_MPa": round(tau_f, 3),
            "clause": "IRC 22:2015 - 606.3.2 Fatigue Strength of Stud Connectors"
        }


    @staticmethod
    def cl_606_4_1_longitudinal_shear_and_spacing(
        V_kN,          # Vertical shear at section (kN)
        beff_mm,       # Effective slab width (mm)
        xu_mm,         # NA depth from top concrete (mm)
        t_slab_mm,     # slab thickness (mm)
        Es=2.0e5,      # MPa
        Ec=30000,      # MPa (secant modulus of concrete)
        As_mm2=0,      # steel beam area mm2
        ys_mm=0,       # CG distance steel to NA
        Ic_mm4=None,   # Composite moment of inertia if already known
        Qu_per_stud_kN=100,  # Shear capacity of 1 stud (kN)
        studs_per_section=2  # default assumption
    ):
        """
        IRC 22 - 606.4.1 Ultimate Limit State Shear Connector Spacing

        Returns:
            VL_kN_per_mm
            spacing_mm
        """

        import math

        # Convert units
        V = V_kN * 1e3   # kN → N
        Qu = Qu_per_stud_kN * 1e3

        # Modular ratio
        n = Es / Ec

        # transformed concrete thickness contributing
        t_eff = min(xu_mm, t_slab_mm)

        # Transformed compressive concrete area
        Aec = n * beff_mm * t_eff  # mm2

        # Y = CG of concrete block from NA
        # xu measured downward.
        # distance from NA to centroid of concrete compressive block
        Y = xu_mm - t_eff / 2

        # Composite inertia if not supplied

        if Ic_mm4 is None:
            Ic_mm4 = (
                (As_mm2 * (ys_mm - xu_mm) ** 2)
                + n * beff_mm * (t_eff ** 3) / 12
                + Aec * (Y ** 2)
            )

        # Longitudinal shear per unit length
        # V_L = V * (Aec * Y) / I
        VL = V * (Aec * Y) / Ic_mm4   # N/mm

        # Spacing
        # Σ Qu / VL
        total_Qu = studs_per_section * Qu

        spacing_mm = total_Qu / VL

        return {
            "modular_ratio_n": round(n, 3),
            "Aec_mm2": round(Aec, 2),
            "Y_mm": round(Y, 2),
            "Ic_mm4": round(Ic_mm4, 2),
            "VL_N_per_mm": round(VL, 3),
            "studs_per_section": studs_per_section,
            "stud_capacity_per_section_kN": round(total_Qu / 1000, 2),
            "spacing_mm": round(spacing_mm, 2),
            "clause": "IRC 22:2015 - 606.4.1 ULS Longitudinal Shear"
        }


    @staticmethod
    def cl_606_4_1_1_full_shear_spacing(
        As_mm2,           # tensile steel area (mm2)
        fy_MPa,           # steel yield strength MPa
        fck_MPa,          # concrete strength MPa
        beff_mm,          # effective slab width
        xu_mm,            # NA depth
        t_slab_mm,        # slab thickness
        Qu_per_stud_kN,   # shear capacity per stud (kN)
        shear_span_mm,    # length between max & zero moment section
        studs_per_section=2,
        gamma_m=1.0
    ):
        """
        IRC 22:2015 Clause 606.4.1.1 Full Shear Connection
        Returns H1, H2, governing H and stud spacing S
        """

        # Effective concrete compressive zone
        t_eff = min(xu_mm, t_slab_mm)
        Aec = beff_mm * t_eff   # mm2

        # Longitudinal forces (kN)
        H1 = (As_mm2 * fy_MPa / gamma_m) * 1e-3
        H2 = (0.36 * fck_MPa * Aec) * 1e-3

        H = min(H1, H2)

        # Total stud capacity in section
        total_Qu_kN = studs_per_section * Qu_per_stud_kN

        # Shear per unit length along span
        H_per_length = H / (shear_span_mm)   # kN/mm

        # Spacing
        spacing_mm = total_Qu_kN / H_per_length

        return {
            "Aec_mm2": round(Aec, 2),
            "H1_kN": round(H1, 3),
            "H2_kN": round(H2, 3),
            "H_governing_kN": round(H, 3),
            "H_per_mm_kN": round(H_per_length, 6),
            "studs_per_section": studs_per_section,
            "total_stud_capacity_kN": total_Qu_kN,
            "spacing_mm": round(spacing_mm, 2),
            "clause": "IRC 22:2015 - 606.4.1.1 Full Shear Connection"
        }


    @staticmethod
    def cl_606_4_2_fatigue_shear_spacing(
        Vr_kN,             # shear range due to LL + impact (kN)
        beff_mm,           
        xu_mm,
        t_slab_mm,
        I_composite_mm4,   # composite moment of inertia
        Qu_per_stud_kN,
        studs_per_section=2
    ):
        """
        IRC 22:2015 Clause 606.4.2
        Serviceability Limit State (Fatigue) Stud Spacing
        """

        # effective concrete depth
        t_eff = min(xu_mm, t_slab_mm)

        # transformed concrete compression area
        Aec = beff_mm * t_eff  # mm2

        # CG distance of slab compression block from NA
        Y = xu_mm - t_slab_mm / 2

        # Longitudinal shear per unit length
        Vr_per_mm = (Vr_kN * Aec * Y) / I_composite_mm4   # kN/mm

        # total stud resistance per section
        total_Qu = studs_per_section * Qu_per_stud_kN

        # spacing
        spacing_mm = total_Qu / Vr_per_mm

        return {
            "Aec_mm2": round(Aec, 2),
            "Y_mm": round(Y, 2),
            "Vr_per_mm_kN": round(Vr_per_mm, 6),
            "total_stud_capacity_kN": total_Qu,
            "spacing_SR_mm": round(spacing_mm, 2),
            "clause": "IRC 22:2015 - 606.4.2 Fatigue Limit State"
        }


    @staticmethod
    def cl_606_6_shear_connector_detailing(
        d_stud_mm,          # stud diameter
        h_stud_mm,          # total stud height
        t_flange_mm,        # thickness of top flange
        edge_distance_mm,   # clear edge distance estud
        ccbottom_mm,        # bottom slab cover to reinforcement
        d_bar_mm,           # bar diameter
        clear_cover_stud_mm=25  # required clear cover to stud top
    ):
        """
        IRC 22:2015 Clause 606.6 - Detailing Requirements Check for Shear Studs
        """

        results = {}

        # 1. Stud diameter limit
        limit_d = 2 * t_flange_mm
        results["stud_diameter_check"] = d_stud_mm <= limit_d
        results["stud_diameter_limit_mm"] = limit_d

        # 2. Minimum stud height
        min_h1 = 4 * d_stud_mm
        min_h2 = 100
        min_required_height = max(min_h1, min_h2)

        results["stud_height_check"] = h_stud_mm >= min_required_height
        results["required_min_height_mm"] = min_required_height

        # 3. Stud head diameter requirement
        min_head_d = 1.5 * d_stud_mm
        results["required_head_diameter_mm"] = min_head_d

        # 4. Edge distance requirement
        results["edge_distance_check"] = edge_distance_mm >= 25
        results["required_edge_distance_mm"] = 25

        # 5. Projection above bottom reinforcement
        min_projection = ccbottom_mm + d_bar_mm + 40
        results["projection_check"] = h_stud_mm >= min_projection
        results["required_min_projection_mm"] = min_projection

        # 6. Overall height with 25 mm cover
        results["overall_height_check"] = h_stud_mm >= 100
        results["required_overall_height_mm"] = 100
        results["required_clear_cover_mm"] = clear_cover_stud_mm

        # Overall status
        results["all_requirements_satisfied"] = all([
            results["stud_diameter_check"],
            results["stud_height_check"],
            results["edge_distance_check"],
            results["projection_check"],
            results["overall_height_check"]
        ])

        results["clause"] = "IRC 22:2015 - Clause 606.6 Detailing of Shear Connectors"

        return results


    @staticmethod
    def cl_606_9_shear_connector_spacing_limits(
        tslab_mm,      # slab thickness
        h_stud_mm,     # stud height
        provided_spacing_mm=None
    ):
        """
        IRC 22:2015 Clause 606.9
        Limiting Criteria for Spacing of Shear Connectors
        """

        # Maximum Allowable Spacing 
        max1 = 600                          # absolute max
        max2 = 3 * tslab_mm                 # 3 * slab thickness
        max3 = 4 * h_stud_mm               # 4 * stud height

        governing_max = min(max1, max2, max3)

        # Minimum Allowable Spacing 
        min_spacing = 75

        result = {
            "max_spacing_limit_mm": governing_max,
            "limit_600_mm": max1,
            "limit_3_tslab_mm": max2,
            "limit_4_hstud_mm": max3,
            "minimum_spacing_limit_mm": min_spacing,
            "clause": "IRC 22:2015 - Clause 606.9 Limiting Criteria for Shear Connector Spacing"
        }

        if provided_spacing_mm is not None:
            result["provided_spacing_mm"] = provided_spacing_mm
            result["is_spacing_acceptable"] = (
                provided_spacing_mm >= min_spacing and
                provided_spacing_mm <= governing_max
            )

        return result

    @staticmethod
    def cl_606_10_transverse_shear_check(
        VL_kN,          # Longitudinal shear per unit length (kN) from earlier clause
        fck,            # MPa
        fyk,            # MPa (yield strength of transverse reinforcement)
        L_mm,           # length of possible shear plane (mm)
        Ast_cm2_per_m,  # provided transverse steel area per metre (cm2/m)
        n_layers=6      # default based on 200mm spacing assumption
    ):
        """
        IRC 22:2015 - Clause 606.10 Transverse Shear Check

        Returns:
            governing_capacity_kN/m
            check_pass
            details
        """

        import math

        # Convert units where required
        VL = VL_kN  # already kN/m
        L = L_mm / 1000.0  # convert mm to metres
        Ast = Ast_cm2_per_m  # cm2/m as required by clause

        # Capacity 1 
        Vcap1 = 0.632 * L * math.sqrt(fck)

        # Capacity 2 
        Vcap2 = 0.232 * L * math.sqrt(fck) + 0.1 * Ast * fyk * n_layers * 1e-3
        # Note: 0.1 * Ast * fyk gives kN since Ast in cm2 and fyk MPa

        # Minimum Steel Requirement 
        Ast_min = (2.5 * VL) / fyk  # cm2/m

        check1 = VL <= Vcap1
        check2 = VL <= Vcap2
        min_steel_ok = Ast >= Ast_min

        return {
            "VL_kN_per_m": round(VL, 3),
            "Vcap1_kN_per_m": round(Vcap1, 3),
            "Vcap2_kN_per_m": round(Vcap2, 3),
            "governing_capacity_kN_per_m": round(max(Vcap1, Vcap2), 3),
            "check_ok": check1 or check2,
            "min_Ast_required_cm2_per_m": round(Ast_min, 3),
            "Ast_provided_ok": min_steel_ok,
            "clause": "IRC 22:2015 - Clause 606.10 Transverse Shear Check"
        }
