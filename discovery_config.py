# Locked Discovery Configuration File
# This defines the exact datasets to exclude from Discovery (held out for Validation)

EXCLUDE_FROM_DISCOVERY = {
    'Kidney': {'GSE30529'},                 # Held out for Validation 2
    'Liver':  {'GSE14323', 'GSE162694'},    # GSE14323 is Val2, GSE162694 is Val1 (sign-inverted)
    'Lungs':  {'GSE83717', 'GSE24206'},     # GSE83717 is Val2, GSE24206 is legacy non-IPF
    'Skin':   {'GSE125362'},                # Held out for Validation 2
}

# Locked Discovery Cohorts per Organ
LOCKED_DISCOVERY_COHORTS = {
    'Kidney': ['GSE104066', 'GSE66494', 'GSE104948', 'GSE104954', 'GSE200818'],
    'Liver':  ['GSE77627', 'GSE89377', 'GSE164760'],
    'Lungs':  ['GSE110147', 'GSE32537', 'GSE53845', 'GSE10667'],
    'Skin':   ['GSE130955', 'GSE95065', 'GSE58095', 'GSE181549'],
}
