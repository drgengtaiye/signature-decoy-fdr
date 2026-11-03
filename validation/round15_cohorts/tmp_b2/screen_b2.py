"""Round 15, part B2: eligibility of the GEO candidates of candidates_r15.csv.

Counts patients and events of tumors of the type with a survival time and an event, from the GEO sample annotations
(external/raw_r15/<GSE>/<GSE>_samples_brief.txt, read as text), with the conventions of
signature-decoy-fdr/external/parse_soft.py: characteristic keys split at the first ':', stripped, duplicates '#2';
filters on lower-cased stripped values (list = membership, string = prefix); time = first number in the value times
time_scale, kept if > 0; event mapped from the lower-cased value. Writes eligibility_r15.csv and cohorts_r15.json.
Run: python3 -I screen_b2.py
"""
import sys, os, re, csv, json, datetime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from soft_meta import load

OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = 1 / 30.4375  # days -> months
W = 7 / 30.4375  # weeks -> months
Y = 12.0         # years -> months
GPL96 = "Affymetrix HG-U133A (GPL96): ~13,000 genes after mapping probes to symbols, below 15,000"

# composite / derived time readers --------------------------------------------------------------------------------
def composite(d, field):
    """'Location: Left; DukesStage: B; ...' stored under one key -> dict of subfields."""
    v = f"{field}: " + d.get(field, "")
    out = {}
    for p in v.split(";"):
        if ":" in p:
            k, x = p.split(":", 1); out[k.strip()] = x.strip()
    return out


def first_num(s):
    m = re.search(r"([-\d.]+)", str(s))
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None


def date_diff_days(a, b):
    try:
        da = datetime.date.fromisoformat(a.strip()); db = datetime.date.fromisoformat(b.strip())
    except Exception:
        return None
    return (db - da).days


def count(spec):
    m = load(spec["gse"])
    n = ev = zero = 0
    for gsm, d in m.items():
        if spec.get("composite"):
            d = dict(d); d.update(composite(d, spec["composite"]))
        ok = True
        for f, allowed in spec.get("filters", {}).items():
            v = str(d.get(f, "nan")).strip().lower()
            ok &= (v in [a.lower() for a in allowed]) if isinstance(allowed, list) else v.startswith(allowed.lower())
        for f, suffix in spec.get("exclude_suffix", {}).items():
            ok &= not str(d.get(f, "")).strip().lower().endswith(suffix.lower())
        if not ok:
            continue
        sets = spec.get("alts") or [(spec["time"], spec["event"], spec["event_map"], spec.get("time_scale", 1.0))]
        for tf, ef, emap, sc in sets:
            if spec.get("time_start"):
                t = date_diff_days(d.get(spec["time_start"], ""), d.get(tf, ""))
                t = None if t is None else t * sc
            else:
                t = first_num(d[tf]) if tf in d else None
                t = None if t is None else t * sc
            e = {str(k).lower(): v for k, v in emap.items()}.get(str(d.get(ef, "")).strip().lower())
            if t is not None and e is not None:
                if t > 0:
                    n += 1; ev += e
                else:
                    zero += 1
                break
    return n, ev, zero


def fmt_map(emap):
    return "; ".join(f"{k}={v}" for k, v in emap.items())


# specs ------------------------------------------------------------------------------------------------------------
# keys: type gse platforms gpl tumor_ok icb endpoint time event event_map time_scale time_unit filters sample_filter
#       platform_flag fail (list of reasons not from counts) notes; optional composite, time_start, exclude_suffix, alts
S = []
def add(**k):
    k.setdefault("tumor_ok", True); k.setdefault("icb", "FALSE"); k.setdefault("filters", {})
    k.setdefault("platform_flag", ""); k.setdefault("fail", []); k.setdefault("notes", "")
    k.setdefault("time_scale", 1.0); k.setdefault("sample_filter", "all samples (all tumors of the type)")
    S.append(k)

NOSURV = "no survival time or event in the GEO sample annotations"

# STAD
add(type="STAD", gse="GSE84437", platforms="GPL6947", gpl="GPL6947", endpoint="OS", time="duration overall survival",
    event="death", event_map={"1": 1, "0": 0}, time_unit="months",
    notes="Round-11 STAD cohort. SuperSeries of GSE84426 (76) and GSE84433 (357), which together are its 433 patients "
          "with survival (all 433 GSMs shared); the other 50 samples (GSE147163) have no survival.")
add(type="STAD", gse="GSE62254", platforms="GPL570", gpl="GPL570", endpoint="none", time="", event="", event_map={},
    time_unit="", fail=[NOSURV], sample_filter="tissue == Gastric tumor (300)",
    notes="ACRG; annotations hold only tissue and patient id (survival is in the paper's supplement, not in GEO). "
          "SubSeries of GSE66229 (same 300 tumors).")
add(type="STAD", gse="GSE15459", platforms="GPL570", gpl="GPL570", endpoint="none", time="", event="", event_map={},
    time_unit="", fail=[NOSURV], sample_filter="tissue == gastric tumor (200)",
    notes="Singapore cohort; annotations hold only tissue (survival in the paper's supplement). SubSeries of GSE15460 "
          "(all 200 GSMs in it).")
add(type="STAD", gse="GSE66229", platforms="GPL570", gpl="GPL570", endpoint="none", time="", event="", event_map={},
    time_unit="", fail=[NOSURV], sample_filter="tissue == Gastric tumor (300 of 400; 100 normal)",
    notes="ACRG SuperSeries of GSE62254 (its 300 tumors) and GSE66222 (100 normals); no survival in GEO.")
add(type="STAD", gse="GSE84433", platforms="GPL6947", gpl="GPL6947", endpoint="OS", time="duration overall survival",
    event="death", event_map={"1": 1, "0": 0}, time_unit="months",
    notes="SubSeries of GSE84437 (round-11 STAD cohort): its 357 patients are among GSE84437's 433 (same patients).")
add(type="STAD", gse="GSE15460", platforms="GPL570; GPL96", gpl=["GPL570", "GPL96"], endpoint="none", time="",
    event="", event_map={}, time_unit="", fail=[NOSURV],
    sample_filter="tissue in {gastric tumor, gastric adenocarcinoma} (287); 70 cell lines and 3 fibroblasts excluded",
    platform_flag="31 samples on GPL96 (U133A, ~13,000 genes)",
    notes="SuperSeries (GSE15455, GSE15456, GSE15459, GSE15537, GSE22183, GSE34942) mixing tumors, cell lines and two "
          "platforms; no survival in GEO annotations. Contains GSE15459.")
add(type="STAD", gse="GSE84426", platforms="GPL6947", gpl="GPL6947", endpoint="OS", time="duration overall survival",
    event="death", event_map={"1": 1, "0": 0}, time_unit="months",
    notes="SubSeries of GSE84437 (round-11 STAD cohort): its 76 patients are among GSE84437's 433.")

# SKCM (mostly metastatic: all melanoma tumors kept)
add(type="SKCM", gse="GSE65904", platforms="GPL10558", gpl="GPL10558", endpoint="DSS",
    time="disease specific survival in days", event="disease specific survival (1=death, 0=alive)",
    event_map={"1": 1, "0": 0}, time_scale=D, time_unit="days",
    sample_filter="all 214 melanoma tumors (metastases and 16 primaries; SKCM is mostly metastatic)",
    notes="Round-11 SKCM cohort. No overall survival annotated; DSS used (DMFS also annotated, 150 with data). "
          "Population-based Lund cohort; no checkpoint treatment annotated.")
add(type="SKCM", gse="GSE54467", platforms="GPL6884", gpl="GPL6884", endpoint="OS",
    time="survival from stage iii tumour banked (months)", event="patient last status",
    event_map={"dead melanoma": 1, "dead not melanoma": 1, "dead cause unknown": 1, "alive nsr": 0,
               "alive with melanoma": 0, "alive status unknown": 0}, time_unit="months",
    sample_filter="all 79 stage III lymph-node melanoma metastases",
    notes="Melanoma Institute Australia stage III nodes, banked before checkpoint therapy era.")
add(type="SKCM", gse="GSE19234", platforms="GPL570", gpl="GPL570", endpoint="OS", time="days lived since metastasis",
    event="staus dead or alive", event_map={"1": 1, "0": 0}, time_scale=D, time_unit="days",
    sample_filter="all 44 metastatic melanoma lesions (some patients contribute several recurrences)",
    notes="44 samples, fewer patients ('B', 'D' suffixes mark repeat lesions of one patient).")
add(type="SKCM", gse="GSE22153", platforms="GPL6102", gpl="GPL6102", endpoint="OS", time="os (days)",
    event="event (0=alive, 1=dead)", event_map={"1": 1, "0": 0}, time_scale=D, time_unit="days",
    sample_filter="all 57 stage IV melanoma metastases", notes="SubSeries of GSE22155 (Lund stage IV, test set).")
add(type="SKCM", gse="GSE59455", platforms="GPL8432", gpl="GPL8432", endpoint="none", time="overall survival (weeks)",
    event="", event_map={}, time_scale=W, time_unit="weeks", fail=["no event (vital status) field in the annotations"],
    sample_filter="all 141 FFPE melanomas (39 primary, 102 metastatic)",
    notes="Annotations give 'overall survival (weeks)' (122 samples) and 'age at death' (134) but no alive/dead status, "
          "so the event cannot be read.")

# LUAD
add(type="LUAD", gse="GSE31210", platforms="GPL570", gpl="GPL570", endpoint="OS", time="days before death/censor",
    event="death", event_map={"dead": 1, "alive": 0}, time_scale=D, time_unit="days",
    filters={"tissue": ["primary lung tumor"]}, sample_filter="tissue == primary lung tumor (226; 20 normal lung out)",
    notes="Stage I-II adenocarcinomas (NCC Tokyo). The authors' flag 'exclude for prognosis analysis due to incomplete "
          "resection or adjuvant therapy' (22 'exclude') is not applied (not part of the rule); with it 204 patients. "
          "Relapse also annotated; with the authors' exclusion 204 patients, 30 deaths (eligible either way).")
add(type="LUAD", gse="GSE72094", platforms="GPL15048", gpl="GPL15048", endpoint="OS", time="survival_time_in_days",
    event="vital_status", event_map={"dead": 1, "alive": 0}, time_scale=D, time_unit="days",
    sample_filter="all 442 lung adenocarcinomas", notes="Round-11 LUAD cohort (Moffitt).")
add(type="LUAD", gse="GSE30219", platforms="GPL570", gpl="GPL570", endpoint="OS", time="follow-up time (months)",
    event="status", event_map={"dead": 1, "alive": 0}, time_unit="months", filters={"histology": ["ADC"]},
    sample_filter="histology == ADC (85 of 293 lung tumors)",
    notes="CIT lung series of mixed histology; LUAD = ADC. Also candidate for LUSC (SQC).")
add(type="LUAD", gse="GSE68465", platforms="GPL96", gpl="GPL96", endpoint="OS", time="months_to_last_contact_or_death",
    event="vital_status", event_map={"dead": 1, "alive": 0}, time_unit="months",
    filters={"disease_state": ["Lung Adenocarcinoma"]}, platform_flag=GPL96, fail=["platform below 15,000 genes (GPL96)"],
    sample_filter="disease_state == Lung Adenocarcinoma (443; 19 normal out)",
    notes="Director's Challenge (Shedden 2008). Would meet patients/events; fails on platform.")
add(type="LUAD", gse="GSE50081", platforms="GPL570", gpl="GPL570", endpoint="OS", time="survival time",
    event="status", event_map={"dead": 1, "alive": 0}, time_scale=Y, time_unit="years",
    filters={"histology": ["adenocarcinoma"]},
    sample_filter="histology == adenocarcinoma (127 of 181 NSCLC; 'NSCLC-favor adenocarcinoma' and adenosquamous out)",
    notes="UHN stage I-II NSCLC.")
add(type="LUAD", gse="GSE37745", platforms="GPL570", gpl="GPL570", endpoint="OS",
    time="days to determined death status", event="dead", event_map={"yes": 1, "no": 0}, time_scale=D,
    time_unit="days", filters={"histology": ["adeno"]}, sample_filter="histology == adeno (106 of 196 NSCLC)",
    notes="Uppsala NSCLC. Also candidate for LUSC (squamous).")
add(type="LUAD", gse="GSE13213", platforms="GPL6480", gpl="GPL6480", endpoint="OS", time="Survival (days)",
    event="Status", event_map={"dead": 1, "alive": 0}, time_scale=D, time_unit="days",
    sample_filter="all 117 lung adenocarcinomas (Histology == AD)",
    notes="Aichi Cancer Center. Shares patients with GSE11969: 87 of GSE11969's 90 adenocarcinomas match a GSE13213 "
          "patient on age, sex and survival days.")
add(type="LUAD", gse="GSE41271", platforms="GPL6884", gpl="GPL6884", endpoint="OS", time="last follow-up survival",
    time_start="date of surgery", event="vital statistics", event_map={"d": 1, "a": 0}, time_scale=D,
    time_unit="days (derived: last follow-up survival - date of surgery)", filters={"histology": ["Adenocarcinoma"]},
    sample_filter="histology == Adenocarcinoma (183 of 275)",
    notes="UNCERTAIN: no survival-time field; time is the difference of two annotated dates (date of surgery, last "
          "follow-up survival), and needs a date-difference reader (parse_soft.py takes a numeric field). MD Anderson; "
          "all 176 samples of GSE42127 are in it (same sample ids). Also candidate for LUSC.")
add(type="LUAD", gse="GSE42127", platforms="GPL6884", gpl="GPL6884", endpoint="OS", time="overall survival months",
    event="survival status", event_map={"d": 1, "a": 0}, time_unit="months",
    filters={"histology": ["Adenocarcinoma"]}, sample_filter="histology == Adenocarcinoma (133 of 176)",
    notes="MD Anderson; all 176 samples (by sample id) are also in GSE41271, so its patients are a subset of GSE41271's.")
add(type="LUAD", gse="GSE26939", platforms="GPL9053", gpl="GPL9053", endpoint="OS", time="survival_months",
    event="survival_status(0='alive',1='dead')", event_map={"1": 1, "0": 0}, time_unit="months",
    sample_filter="all 116 lung adenocarcinomas ('squamoid' is an expression subtype of adenocarcinoma)",
    notes="UNC two-colour Agilent custom 4x44K (GPL9053, 33,421 features); symbol annotation to be checked after "
          "download. SubSeries of GSE36471.")
add(type="LUAD", gse="GSE3141", platforms="GPL570", gpl="GPL570", endpoint="OS", time="Surv(months)",
    event="STATUS(0=alive, 1=dead)", event_map={"1": 1, "0": 0}, time_unit="months",
    composite_per_title=True, filters={"Cell type": ["A"]},
    sample_filter="Cell type == A (adenocarcinoma) of 111 NSCLC",
    notes="Annotations are one line per sample keyed by the sample title ('Cell type: A; Surv(months): ...; STATUS...').")
add(type="LUAD", gse="GSE11969", platforms="GPL7015", gpl="GPL7015", endpoint="OS", time="Survival (days)",
    event="Status", event_map={"dead": 1, "alive": 0}, time_scale=D, time_unit="days", filters={"Histology": ["AD"]},
    sample_filter="Histology == AD (90 of 163)",
    platform_flag="custom Agilent 21.6K array (21,619 features; authors report 18,175 unique genes); may fall near "
                  "15,000 symbols, check after download",
    notes="Aichi Cancer Center; 87 of its 90 adenocarcinomas match GSE13213 patients on age, sex and survival days "
          "(same patients).")
add(type="LUAD", gse="GSE29016", platforms="GPL6947", gpl="GPL6947", endpoint="OS", time="os", event="osbin",
    event_map={"1": 1, "0": 0}, time_unit="years (presumed)", filters={"histology": ["AC"]},
    sample_filter="histology == AC (38 of 72)", notes="Lund early-stage lung carcinomas of mixed histology.")
add(type="LUAD", gse="GSE29013", platforms="GPL570", gpl="GPL570", endpoint="OS", time="death_time",
    event="death_event", event_map={"1": 1, "0": 0}, time_scale=Y, time_unit="years",
    filters={"histology": ["Adenocarcinoma"]}, sample_filter="histology == Adenocarcinoma (30 of 55)",
    notes="MD Anderson FFPE NSCLC.")

# LUSC
add(type="LUSC", gse="GSE73403", platforms="GPL6480", gpl="GPL6480", endpoint="OS", time="survival (year)",
    event="events (death = 1, alive = 0)", event_map={"1": 1, "0": 0}, time_scale=Y, time_unit="years",
    sample_filter="all 69 lung squamous carcinomas")
add(type="LUSC", gse="GSE30219", platforms="GPL570", gpl="GPL570", endpoint="OS", time="follow-up time (months)",
    event="status", event_map={"dead": 1, "alive": 0}, time_unit="months", filters={"histology": ["SQC"]},
    sample_filter="histology == SQC (61 of 293 lung tumors)", notes="CIT lung series; also candidate for LUAD.")
add(type="LUSC", gse="GSE74777", platforms="GPL17586", gpl="GPL17586", endpoint="OS", time="days before death/censor",
    event="death", event_map={"dead": 1, "alive": 0}, time_scale=D, time_unit="days",
    sample_filter="all 107 lung squamous cell carcinomas", notes="NCC Tokyo stage I-II, Affymetrix HTA 2.0.")
add(type="LUSC", gse="GSE37745", platforms="GPL570", gpl="GPL570", endpoint="OS",
    time="days to determined death status", event="dead", event_map={"yes": 1, "no": 0}, time_scale=D,
    time_unit="days", filters={"histology": ["squamous"]}, sample_filter="histology == squamous (66 of 196)",
    notes="Uppsala NSCLC; also candidate for LUAD.")
add(type="LUSC", gse="GSE41271", platforms="GPL6884", gpl="GPL6884", endpoint="OS", time="last follow-up survival",
    time_start="date of surgery", event="vital statistics", event_map={"d": 1, "a": 0}, time_scale=D,
    time_unit="days (derived: last follow-up survival - date of surgery)", filters={"histology": ["Squamous"]},
    sample_filter="histology == Squamous (80 of 275; Sarcomatoid-squamous out)",
    notes="UNCERTAIN: time is the difference of two annotated dates (see LUAD row). Also candidate for LUAD.")
add(type="LUSC", gse="GSE157011", platforms="GPL570", gpl="GPL570", endpoint="OS", time="os_mo", event="os_event",
    event_map={"1": 1, "0": 0}, time_unit="months", sample_filter="all 484 squamous cell carcinomas of lung",
    notes="Round-11 LUSC cohort (Alliance). SuperSeries of GSE157009 and GSE157010.")

# COAD (colon only: rectum and rectosigmoid are TCGA READ)
COLON_103479 = ["Sigmoid", "Right colon", "Left Colon", "Caecum", "Ascending colon", "Hepatic flexure",
                "Transverse colon", "Asc colon", "colon", "Left colon", "Colon", "right colon",
                "Transverse/left colon", "left colon/cecum", "sigma", "Descending colon", "Caecum + Asc colon"]
add(type="COAD", gse="GSE39582", platforms="GPL570", gpl="GPL570", endpoint="OS", time="os.delay (months)",
    event="os.event", event_map={"1": 1, "0": 0}, time_unit="months",
    filters={"dataset": ["discovery", "validation"]}, sample_filter="dataset in {discovery, validation} (566 tumors; 19 "
    "non-tumoral out)", notes="Round-11 COAD cohort (CIT colon cancer). SubSeries of GSE40967.")
add(type="COAD", gse="GSE17536", platforms="GPL570", gpl="GPL570", endpoint="OS",
    time="overall survival follow-up time", event="overall_event (death from any cause)",
    event_map={"death": 1, "no death": 0}, time_unit="months", sample_filter="all 177 colon cancers (Moffitt)",
    notes="SubSeries of GSE17538 (same 177 patients). GSE29621 maps 63 of its 65 samples to GSE17536 samples, and "
          "GSE14333 shares those patients (same tumor ids).")
add(type="COAD", gse="GSE17538", platforms="GPL570; GPL1261", gpl="GPL570", endpoint="OS",
    time="overall survival follow-up time", event="overall_event (death from any cause)",
    event_map={"death": 1, "no death": 0}, time_unit="months",
    filters={"(source)": ["Moffitt Cancer Center", "Vanderbilt University Medical Center"]},
    sample_filter="(source) in {Moffitt Cancer Center, Vanderbilt University Medical Center} (232 colon cancers); "
                  "6 adenomas and 6 mouse MC38 samples (GPL1261) out",
    notes="SuperSeries = GSE17536 (177) + GSE17537 (55): same patients as those two (all 232 GSMs shared). The mouse "
          "GPL1261 samples are a separate experiment, not a cohort.")
add(type="COAD", gse="GSE38832", platforms="GPL570", gpl="GPL570", endpoint="DSS",
    time="dss_time (disease specific survival time, months)", event="dss_event (disease specific survival)",
    event_map={"1 (death from cancer)": 1, "0 (no death)": 0}, time_unit="months",
    sample_filter="all 122 colon cancers (Vanderbilt 104, UAB 14, Nashville VA 4)",
    notes="No overall survival annotated; DSS used (DFS for 92). Vanderbilt patients may overlap GSE17537 "
          "(same group); not checkable from the annotations.")
add(type="COAD", gse="GSE17537", platforms="GPL570", gpl="GPL570", endpoint="OS",
    time="overall survival follow-up time", event="overall_event (death from any cause)",
    event_map={"death": 1, "no death": 0}, time_unit="months", sample_filter="all 55 colon cancers (Vanderbilt)",
    notes="SubSeries of GSE17538.")
add(type="COAD", gse="GSE29621", platforms="GPL570", gpl="GPL570", endpoint="DFS",
    time="disease-free survival time (dfs)", event="dfs event", event_map={"recurrence": 1, "no recurrence": 0},
    time_unit="months", sample_filter="all 65 colon cancers",
    notes="OS time annotated without an OS event; DFS only. 63 of 65 samples map to GSE17536 samples (field "
          "'corresponding to gse17536'): same patients.")
add(type="COAD", gse="GSE103479", platforms="GPL23985", gpl="GPL23985", endpoint="OS", time="overall survival time",
    event="status alive.dead", event_map={"dead": 1, "alive": 0}, time_unit="months",
    filters={"tumour site": COLON_103479},
    sample_filter="tumour site is colon (128 of 156; rectum, upper rectum, rectum intraperitoneal, rectosigmoid, "
                  "recto-sigmoid and 'Distal' (ambiguous) out)",
    notes="Stage II-III CRC from 4 European centres; GPL23985 is a custom array (110,961 features), symbol "
          "annotation to be checked after download. Adjuvant 5FU in 66 (not checkpoint). All sites: 155 patients, 60 deaths.")
add(type="COAD", gse="GSE14333", platforms="GPL570", gpl="GPL570", endpoint="DFS", time="DFS_Time", event="DFS_Cens",
    event_map={"0": 1, "1": 0}, time_unit="months", composite="Location",
    filters={"Location": ["Left", "Right", "Colon"]},
    sample_filter="Location in {Left, Right, Colon} (250 of 290; 39 Rectum and 1 blank out); Dukes D have no DFS",
    notes="UNCERTAIN coding: all annotations sit in one line ('Location: Left; DukesStage: B; ...; DFS_Time: 22.38; "
          "DFS_Cens: 0; ...'), needing a composite reader. DFS_Cens read as 1 = censored, 0 = recurrence (median DFS "
          "19.1 months for 0 vs 47.6 for 1; 1 = event would give 85% recurrence in Dukes B). Only DFS annotated. "
          "Shares patients with GSE17536 (63 tumor ids also in GSE29621, which maps them to GSE17536). All sites incl. rectum: 226/50.")
add(type="COAD", gse="GSE72970", platforms="GPL570", gpl="GPL570", endpoint="OS", time="os", event="os censored",
    event_map={"1": 1, "0": 0}, time_unit="months",
    filters={"tumor location": ["Left colon", "Right colon", "Caecum", "Transverse colon"]},
    sample_filter="tumor location is colon (86 of 124; 19 rectum and 19 rectum-sigmoid junction out)",
    notes="UNCERTAIN coding: 'os censored' read as 1 = death, 0 = censored (median OS 18.0 months for 1 vs 47.9 for 0; "
          "92/124 deaths fits synchronous metastatic CRC on first-line chemotherapy). Primary tumors of metastatic "
          "patients; chemotherapy with bevacizumab or cetuximab in some (not checkpoint). SuperSeries of GSE72968, "
          "GSE72969. Counting the rectosigmoid junction as colon: 105 patients, 78 deaths; all sites 124/92.")
add(type="COAD", gse="GSE87211", platforms="GPL13497", gpl="GPL13497", endpoint="DSS", time="survival time (month)",
    event="death due to tumor", event_map={"1": 1, "0": 0}, time_unit="months", filters={"tissue": ["rectal tumor"]},
    tumor_ok=False, fail=["rectal cancer (TCGA READ), not colon (COAD)"],
    sample_filter="tissue == rectal tumor (203; 160 mucosa out)",
    notes="Rectal cancers after preoperative radiochemotherapy; only tumor-related death annotated (DSS).")
add(type="COAD", gse="GSE40967", platforms="GPL570; GPL16070", gpl="GPL570", endpoint="OS", time="os.delay (months)",
    event="os.event", event_map={"1": 1, "0": 0}, time_unit="months",
    filters={"(platform)": ["GPL570"], "dataset": ["discovery", "validation"]},
    sample_filter="(platform) == GPL570 and dataset in {discovery, validation} (566 tumors); 463 CGH arrays "
                  "(GPL16070) out",
    notes="SuperSeries of GSE39582 (expression, the round-11 COAD cohort; all 585 GSMs shared) and GSE40966 (CGH "
          "arrays, not expression). Its expression samples ARE GSE39582: same patients, same arrays.")

# CESC
add(type="CESC", gse="GSE44001", platforms="GPL14951", gpl="GPL14951", endpoint="DFS",
    time="disease_free_survival_(dfs)_(months)", event="status_of_dfs", event_map={"1": 1, "0": 0},
    time_unit="months", sample_filter="all 300 early cervical cancers", notes="Round-11 CESC cohort. Only DFS annotated.")
add(type="CESC", gse="GSE52903", platforms="GPL6244", gpl="GPL6244", endpoint="none", time="", event="", event_map={},
    time_unit="", fail=[NOSURV], sample_filter="clinic == Case (55 of 72)",
    notes="Only stage, histology and age annotated. SubSeries of GSE52904.")
add(type="CESC", gse="GSE30759", platforms="GPL8490", gpl="GPL8490", endpoint="OS", time="survival time (years)",
    event="death status", event_map={"1": 1, "0": 0}, time_scale=Y, time_unit="years",
    filters={"status(0=normal,1=cancer)": ["1"]}, fail=["DNA methylation array (Infinium 27k), not expression"],
    platform_flag="GPL8490 is a DNA methylation array: no expression",
    sample_filter="status(0=normal,1=cancer) == 1 (48 of 63)", notes="SubSeries of GSE30760.")
add(type="CESC", gse="GSE39001", platforms="GPL201; GPL6244", gpl=["GPL201", "GPL6244"], endpoint="none", time="",
    event="", event_map={}, time_unit="", fail=[NOSURV], sample_filter="clinic == Case (62 of 79)",
    platform_flag="55 samples on GPL201 (HG-Focus, ~8,700 genes)",
    notes="Two platforms (GPL201 55, GPL6244 24); no survival annotated.")

# SARC
add(type="SARC", gse="GSE21050", platforms="GPL570", gpl="GPL570", endpoint="MFS", time="time", event="metastasis",
    event_map={"yes": 1, "no": 0}, time_scale=Y, time_unit="years",
    sample_filter="all 310 untreated primary soft-tissue sarcomas",
    notes="Round-11 SARC cohort. Only metastasis-free survival annotated. 303 of its sample ids (S#) recur in GSE71118.")
add(type="SARC", gse="GSE17674", platforms="GPL570", gpl="GPL570", endpoint="OS", time="ovs (months)", event="status",
    event_map={"dead": 1, "ned": 0, "awd": 0}, time_unit="months", filters={"state": ["Primary"]}, tumor_ok=False,
    fail=["Ewing sarcoma family, not a TCGA SARC histology"],
    sample_filter="state == Primary (32 of 44 Ewing tumors; 18 normal muscle out)",
    notes="Ewing sarcoma family (paediatric bone/soft tissue); TCGA SARC holds adult soft-tissue sarcomas. "
          "SubSeries of GSE17679.")
add(type="SARC", gse="GSE30929", platforms="GPL96", gpl="GPL96", endpoint="MFS", time="tt.drfs", event="drfs",
    event_map={"true": 1, "false": 0}, time_unit="months", platform_flag=GPL96,
    fail=["platform below 15,000 genes (GPL96)"], sample_filter="all 140 primary liposarcomas",
    notes="Distant recurrence-free survival only (treated as MFS).")
add(type="SARC", gse="GSE63157", platforms="GPL5175", gpl="GPL5175", endpoint="OS", time="survival time (days)",
    event="survival status", event_map={"dead": 1, "alive": 0}, time_scale=D, time_unit="days",
    alts=[("survival time (days)", "survival status", {"dead": 1, "alive": 0}, D),
          ("os-time", "os (0=alive; 1=dead)", {"1": 1, "0": 0}, D)],
    tumor_ok=False, fail=["Ewing sarcoma, not a TCGA SARC histology"],
    sample_filter="all 85 primary Ewing sarcomas (two subseries with different field names: 46 + 39)",
    notes="COG Ewing sarcoma (SuperSeries of GSE63155, GSE63156); survival fields differ between the subseries.")
add(type="SARC", gse="GSE71118", platforms="GPL570", gpl="GPL570", endpoint="MFS", time="time", event="metastasis",
    event_map={"yes": 1, "no": 0}, time_scale=Y, time_unit="years",
    sample_filter="all 312 frozen sarcomas (undifferentiated, leiomyo-, dedifferentiated lipo-, myxofibro-, other)",
    notes="Only MFS annotated. Largely the same patients as GSE21050 (303 of 312 sample ids S# shared; French "
          "Sarcoma Group CINSARC series). SubSeries of GSE71121.")

# HNSC
add(type="HNSC", gse="GSE65858", platforms="GPL10558", gpl="GPL10558", endpoint="OS", time="os", event="os_event",
    event_map={"true": 1, "false": 0}, time_scale=D, time_unit="days", filters={"tumor_type": ["Primary"]},
    sample_filter="tumor_type == Primary (253 of 270; 14 secondary, 3 relapse out)", notes="Round-11 HNSC cohort.")
add(type="HNSC", gse="GSE41613", platforms="GPL570", gpl="GPL570", endpoint="OS", time="fu time", event="vital",
    event_map={"dead-oral ca": 1, "dead-non oc": 1, "dead -unk cause": 1, "alive": 0}, time_unit="months",
    sample_filter="all 97 HPV-negative oral squamous cell carcinomas",
    notes="Any death counted (oral-cancer death 30, other 14, unknown 7).")
add(type="HNSC", gse="GSE27020", platforms="GPL96", gpl="GPL96", endpoint="DFS", time="dfs (months)",
    event="dfs status (1 = recurred)", event_map={"1": 1, "0": 0}, time_unit="months", platform_flag=GPL96,
    fail=["platform below 15,000 genes (GPL96)"], sample_filter="all 109 laryngeal tumors")
add(type="HNSC", gse="GSE42743", platforms="GPL570", gpl="GPL570", endpoint="OS", time="futime", event="alldeath",
    event_map={"1": 1, "0": 0}, time_unit="unclear (days or months)", filters={"tu-vs-nl": ["cancer"]},
    sample_filter="tu-vs-nl == cancer (74 of 103; 29 adjacent normal out)", notes="MD Anderson oral cavity.")

# BLCA
add(type="BLCA", gse="GSE13507", platforms="GPL6102", gpl="GPL6102", endpoint="OS", time="survival month",
    event="overall survival", event_map={"death": 1, "survival": 0}, time_unit="months",
    filters={"biological source": ["primary bladder cancer"]},
    sample_filter="biological source == primary bladder cancer (165 of 256)", notes="Round-11 BLCA cohort.")
add(type="BLCA", gse="GSE31684", platforms="GPL570", gpl="GPL570", endpoint="OS", time="survival.months",
    event="last known status", event_map={"dod": 1, "doc": 1, "ned": 0}, time_unit="months",
    sample_filter="all 93 cystectomy urothelial carcinomas",
    notes="OS event = dead of disease (38) or of other causes (27); NED = alive.")
add(type="BLCA", gse="GSE32894", platforms="GPL6947", gpl="GPL6947", endpoint="DSS", time="time_to_dod_(months)",
    event="dod_event_(yes/no)", event_map={"yes": 1, "no": 0}, time_unit="months",
    sample_filter="all 308 urothelial carcinomas; 224 with disease-specific follow-up (80 'dod_excluded')",
    notes="No overall survival annotated; DSS used. Exactly at the 25-event threshold. GSE32548 (131) is a subset "
          "(111 of its UC ids in GSE32894).")
add(type="BLCA", gse="GSE48075", platforms="GPL6947", gpl="GPL6947", endpoint="OS", time="survival (mo)",
    event="os censor", event_map={"uncensored": 1, "censored": 0}, time_unit="months",
    sample_filter="all 142 primary bladder tumors; survival for the 73 muscle-invasive",
    notes="MD Anderson; SubSeries of GSE48277.")
add(type="BLCA", gse="GSE32548", platforms="GPL6947", gpl="GPL6947", endpoint="none", time="", event="", event_map={},
    time_unit="", fail=[NOSURV], sample_filter="all 131 bladder tumors",
    notes="Only stage and grade annotated; 111 of its tumors are in GSE32894.")
add(type="BLCA", gse="GSE48276", platforms="GPL14951", gpl="GPL14951", endpoint="OS", time="survival (mo)",
    event="os censor", event_map={"uncensored": 1, "censored": 0}, time_unit="months",
    sample_filter="all 116 FFPE bladder tumors; survival for 73", notes="MD Anderson FFPE; SubSeries of GSE48277.")

# OV (TCGA OV = serous cystadenocarcinoma: serous kept where histology is annotated)
add(type="OV", gse="GSE26193", platforms="GPL570", gpl="GPL570", endpoint="OS", time="os time (years)",
    event="os event", event_map={"1": 1, "0": 0}, time_scale=Y, time_unit="years",
    filters={"histological type": ["Serous"]}, sample_filter="histological type == Serous (79 of 107)",
    notes="UNCERTAIN: fails only because non-serous carcinomas (endometrioid 8, mucinous 8, clear cell 6, other 6) are "
          "left out as not TCGA OV; with all 107 ovarian carcinomas it would have 107 patients, 76 deaths.")
add(type="OV", gse="GSE140082", platforms="GPL14951", gpl="GPL14951", endpoint="OS", time="final_ostm",
    event="final_osid", event_map={"1": 1, "0": 0}, time_scale=D, time_unit="days",
    filters={"histology.serous": ["serous"]}, sample_filter="histology.serous == serous (277 of 380; 'other' out)",
    notes="ICON7 trial FFPE (carboplatin-paclitaxel with or without bevacizumab; not checkpoint). With all 380 "
          "carcinomas 380 patients, 96 deaths (eligible either way).")
add(type="OV", gse="GSE26712", platforms="GPL96", gpl="GPL96", endpoint="DSS", time="survival years", event="status",
    event_map={"dod (dead of disease)": 1, "ned (no evidence of disease)": 0, "awd (alive with disease)": 0},
    time_scale=Y, time_unit="years", filters={"tissue": ["Late-stage high-grade ovarian cancer"]},
    platform_flag=GPL96, fail=["platform below 15,000 genes (GPL96)"],
    sample_filter="tissue == Late-stage high-grade ovarian cancer (185; 10 normal OSE out)",
    notes="Only death of disease annotated (DSS).")
add(type="OV", gse="GSE63885", platforms="GPL570", gpl="GPL570", endpoint="OS", time="os - overall survival [days]",
    event="clinical status at last follow-up (awd - alive with disease; dod - dead from disease; ned - no evidence of disease)",
    event_map={"dod": 1, "ned": 0, "awd": 0}, time_scale=D, time_unit="days",
    filters={"histophatological type of tumor": ["serous ovarian ca."]},
    sample_filter="histophatological type of tumor == serous ovarian ca. (73 of 101)",
    notes="Deaths recorded only as 'dead from disease'. With all histologies 75 patients, 66 deaths: still below 80.")
add(type="OV", gse="GSE9891", platforms="GPL570", gpl="GPL570", endpoint="none", time="", event="", event_map={},
    time_unit="", fail=[NOSURV], sample_filter="Type == Malignant, Subtype == Ser/PapSer (would be the filter)",
    notes="AOCS; annotations hold site, type, subtype, stage, grade only (survival in the paper's supplement). "
          "SubSeries of GSE9899.")
add(type="OV", gse="GSE14764", platforms="GPL96", gpl="GPL96", endpoint="OS", time="overall survival time",
    event="overall survival event", event_map={"1": 1, "0": 0}, time_unit="months",
    filters={"histological type": ["serous ovca"]}, platform_flag=GPL96, fail=["platform below 15,000 genes (GPL96)"],
    sample_filter="histological type == serous ovca (68 of 80)", notes="With all 80: 21 deaths, also below 25.")
add(type="OV", gse="GSE30161", platforms="GPL570", gpl="GPL570", endpoint="OS", time="overall survival days",
    event="censoring(dead=1,alive=0)", event_map={"1": 1, "0": 0}, time_scale=D, time_unit="days",
    filters={"histo": ["Serous"]}, sample_filter="histo == Serous (47 of 58)", notes="FFPE, stage III-IV.")
add(type="OV", gse="GSE53963", platforms="GPL6480", gpl="GPL6480", endpoint="none", time="", event="", event_map={},
    time_unit="", fail=[NOSURV], sample_filter="all 174 high-grade serous tumors",
    notes="Only 'tissue: ovarian' annotated; two-colour against a reference pool.")
add(type="OV", gse="GSE32062", platforms="GPL6480; GPL570", gpl="GPL6480", endpoint="OS", time="os (m)",
    event="death (1)", event_map={"1": 1, "0": 0}, time_unit="months", filters={"(platform)": ["GPL6480"]},
    sample_filter="(platform) == GPL6480 (260 high-grade serous; 10 GPL570 cross-platform repeats out)",
    notes="Round-11 OV cohort.")

# BRCA
add(type="BRCA", gse="GSE20685", platforms="GPL570", gpl="GPL570", endpoint="OS", time="follow_up_duration (years)",
    event="event_death", event_map={"1": 1, "0": 0}, time_scale=Y, time_unit="years",
    sample_filter="all 327 primary breast cancers", notes="Taiwan; 30 had neoadjuvant chemotherapy (not checkpoint).")
add(type="BRCA", gse="GSE42568", platforms="GPL570", gpl="GPL570", endpoint="OS", time="overall survival time_days",
    event="overall survival event", event_map={"1": 1, "0": 0}, time_scale=D, time_unit="days",
    filters={"tissue": ["breast cancer"]}, sample_filter="tissue == breast cancer (104 of 121; 17 normal out)",
    notes="Dublin, untreated biopsies.")
add(type="BRCA", gse="GSE58812", platforms="GPL570", gpl="GPL570", endpoint="OS", time="os (days)", event="death",
    event_map={"1": 1, "0": 0}, time_scale=D, time_unit="days", sample_filter="all 107 triple-negative breast cancers",
    notes="SubSeries of GSE103091 (its 107 GSMs are GSE103091's 107 samples with survival: same patients).")
add(type="BRCA", gse="GSE21653", platforms="GPL570", gpl="GPL570", endpoint="DFS", time="dfs time (months)",
    event="dfs evt", event_map={"1": 1, "0": 0}, time_unit="months", sample_filter="all 266 early breast cancers",
    notes="Institut Paoli-Calmettes; only DFS annotated.")
add(type="BRCA", gse="GSE20711", platforms="GPL570", gpl="GPL570", endpoint="OS", time="t.os", event="e.os",
    event_map={"1": 1, "0": 0}, time_scale=Y, time_unit="years (values end in ' y')",
    filters={"(source)": ["Breast tumor"]}, sample_filter="(source) == Breast tumor (88 of 90; 2 normal out)",
    notes="Exactly at the 25-event threshold. SubSeries of GSE20713.")
add(type="BRCA", gse="GSE96058", platforms="GPL11154; GPL18573", gpl=["GPL11154", "GPL18573"], endpoint="OS",
    time="overall survival days", event="overall survival event", event_map={"1": 1, "0": 0}, time_scale=D,
    time_unit="days", exclude_suffix={"(title)": "repl"},
    sample_filter="all primary breast tumors, technical replicates (title ending 'repl', 136) out: 3,273 patients",
    notes="Round-11 BRCA cohort (SCAN-B; RNA-seq on HiSeq 2000 and NextSeq 500, one cohort; expression from the "
          "round-7 preparation of the supplementary matrix, not SOFT tables).")
add(type="BRCA", gse="GSE7390", platforms="GPL96", gpl="GPL96", endpoint="OS", time="t.os", event="e.os",
    event_map={"1": 1, "0": 0}, time_scale=D, time_unit="days", platform_flag=GPL96,
    fail=["platform below 15,000 genes (GPL96)"], sample_filter="all 198 node-negative untreated breast cancers",
    notes="TRANSBIG.")
add(type="BRCA", gse="GSE103091", platforms="GPL570", gpl="GPL570", endpoint="OS", time="os (days)", event="death",
    event_map={"1": 1, "0": 0}, time_scale=D, time_unit="days",
    sample_filter="all 238 TNBC; survival only for the 107 of GSE58812 (GSE83937's 131 have none)",
    notes="SuperSeries of GSE58812 and GSE83937; the patients with survival are exactly GSE58812's 107 (same GSMs).")


def count_any(s):
    if s["endpoint"] == "none" or not s.get("event"):
        if s["gse"] == "GSE59455":
            return 0, 0, 0
        return 0, 0, 0
    if s.get("composite_per_title"):
        # GSE3141: one characteristics line keyed by the title, value 'Cell type: A; Surv(months): ...; STATUS...: 0'
        m = load(s["gse"]); n = ev = zero = 0
        for gsm, d in m.items():
            v = d.get(d["(title)"], "")
            kv = {}
            for p in v.split(";"):
                if ":" in p:
                    k, x = p.split(":", 1); kv[k.strip()] = x.strip()
            if kv.get("Cell type", "").lower() not in [a.lower() for a in s["filters"]["Cell type"]]:
                continue
            t = first_num(kv.get(s["time"], "")); e = {k.lower(): v for k, v in s["event_map"].items()}.get(
                kv.get(s["event"], "").strip().lower())
            if t is not None and e is not None:
                if t > 0:
                    n += 1; ev += e
                else:
                    zero += 1
        return n, ev, zero
    return count(s)


rows, cfg = [], {}
for s in S:
    n, ev, zero = count_any(s)
    reasons = list(s["fail"])
    if not s["tumor_ok"] and not any("TCGA" in r or "rectal" in r for r in reasons):
        reasons.append("not tumors of the type")
    if s["icb"] == "TRUE":
        reasons.append("checkpoint-treated")
    if s["endpoint"] != "none":
        if n < 80:
            reasons.append(f"fewer than 80 patients ({n})")
        if ev < 25:
            reasons.append(f"fewer than 25 events ({ev})")
    elig = not reasons
    notes = s["notes"]
    if zero:
        notes = (notes + " " if notes else "") + f"{zero} sample(s) with time <= 0 not counted (loader keeps time > 0)."
    rows.append({"type": s["type"], "gse": s["gse"], "platforms": s["platforms"],
                 "tumor_ok": "TRUE" if s["tumor_ok"] else "FALSE", "icb_treated": s["icb"], "endpoint": s["endpoint"],
                 "time_field": s["time"] if not s.get("time_start") else f"{s['time']} minus {s['time_start']}",
                 "event_field": s["event"], "event_coding": fmt_map(s["event_map"]), "time_unit": s["time_unit"],
                 "sample_filter": s["sample_filter"], "n_patients": n, "n_events": ev,
                 "platform_flag": s["platform_flag"], "eligible": "TRUE" if elig else "FALSE",
                 "reason": "; ".join(reasons), "notes": notes})
    if elig:
        c = {"gse": s["gse"], "gpl": s["gpl"], "endpoint": s["endpoint"], "time": s["time"], "event": s["event"],
             "event_map": {k.lower(): v for k, v in s["event_map"].items()}, "time_scale": s["time_scale"]}
        if s["filters"]:
            c["filters"] = s["filters"]
        if s.get("time_start"):
            c["time_start"] = s["time_start"]
            c["time_kind"] = "date_difference_days (time = date in 'time' minus date in 'time_start', ISO dates)"
        if s.get("composite"):
            c["composite"] = s["composite"]
            c["composite_note"] = ("all annotations are in one characteristics line keyed '%s'; parse '%s: <value>' "
                                   "split on ';' into key: value pairs; time, event and filters refer to those keys"
                                   % (s["composite"], s["composite"]))
        if s.get("exclude_suffix"):
            c["exclude_suffix"] = s["exclude_suffix"]
        c["n_patients"] = n; c["n_events"] = ev
        cfg[f"{s['type']}_{s['gse']}"] = c

cols = ["type", "gse", "platforms", "tumor_ok", "icb_treated", "endpoint", "time_field", "event_field", "event_coding",
        "time_unit", "sample_filter", "n_patients", "n_events", "platform_flag", "eligible", "reason", "notes"]
with open(os.path.join(OUT, "eligibility_r15.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
with open(os.path.join(OUT, "cohorts_r15.json"), "w") as f:
    json.dump(cfg, f, indent=1)
for r in rows:
    print(f"{r['type']:5s} {r['gse']:10s} {r['endpoint']:4s} n={r['n_patients']:5d} ev={r['n_events']:4d} "
          f"elig={r['eligible']:5s} {r['reason']}")
print("pairs", len(rows), "eligible", sum(r["eligible"] == "TRUE" for r in rows))
