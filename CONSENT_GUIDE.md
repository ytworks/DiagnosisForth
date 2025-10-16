# Consent and Terms of Service Guide

## Overview

DiagnosisForth requires explicit acceptance of the Terms of Service before use. This ensures users are aware of:
- Data collection for improvement and debugging
- Complete exemption of author liability
- Assumption of all risks by the user

## First Time Setup

### In Jupyter Notebook

1. Import the consent module:
```python
from DiagnosisForth.consent import accept_terms
```

2. Run the consent function:
```python
accept_terms()
```

3. Fill in your name and click "I Accept the Terms"

4. Restart the kernel

5. Import DiagnosisForth normally:
```python
import DiagnosisForth
```

### In Terminal/Command Line

When you first import DiagnosisForth, you will be prompted to:
1. Read the terms of service
2. Enter your full name as a digital signature
3. Type 'yes' to accept the terms

## Consent Management Functions

### Check Consent Status
```python
from DiagnosisForth.consent import check_consent_status
check_consent_status()
```

### Reset Consent (for testing)
```python
from DiagnosisForth.consent import reset_consent
reset_consent()
```

## Important Notes

- Consent is stored locally in `~/.diagnosisforth/consent.json`
- Consent is version-specific (you may need to re-accept for major updates)
- Your acceptance is logged for compliance purposes
- IP address and timestamp are recorded with your consent

## Privacy and Data Collection

When you accept the terms:
- Your name, IP address, and acceptance timestamp are logged
- Usage patterns may be monitored for improvement purposes
- Error reports may be collected for debugging

All data collection only occurs AFTER explicit consent is given.

## Troubleshooting

### "Terms not accepted" error in Jupyter
1. Run `accept_terms()` in a new cell
2. Restart the kernel after accepting
3. Import the module again

### Consent form not showing in Jupyter
- Ensure ipywidgets is installed: `pip install ipywidgets`
- If widgets don't work, the text-based consent will be used

### Need to re-accept terms
```python
from DiagnosisForth.consent import reset_consent, accept_terms
reset_consent()
accept_terms()
```