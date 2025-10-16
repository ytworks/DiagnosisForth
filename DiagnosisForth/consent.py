"""
Consent management module for DiagnosisForth.
This module provides functions to manage user consent for terms of service.
"""

from .utils.terms_consent import (
    check_existing_consent,
    request_consent,
    request_consent_jupyter,
    is_jupyter
)


def accept_terms():
    """
    Interactive function to accept terms of service in Jupyter notebooks.
    Call this function in a Jupyter cell to display the consent form.
    
    Example:
        >>> from DiagnosisForth.consent import accept_terms
        >>> accept_terms()
    """
    if not check_existing_consent():
        if is_jupyter():
            request_consent_jupyter()
        else:
            request_consent()
    else:
        print("✅ You have already accepted the terms of service.")
        print("You can use the DiagnosisForth module.")


def reset_consent():
    """
    Reset consent status (for testing or re-consent).
    This will delete the stored consent and require re-acceptance.
    """
    from pathlib import Path
    consent_file = Path.home() / ".diagnosisforth" / "consent.json"
    
    if consent_file.exists():
        consent_file.unlink()
        print("⚠️ Consent has been reset.")
        print("You will need to accept the terms again to use the software.")
    else:
        print("No existing consent found.")


def check_consent_status():
    """
    Check the current consent status.
    Returns information about whether terms have been accepted.
    """
    consent = check_existing_consent()
    
    if consent:
        print("✅ Terms of Service Status: ACCEPTED")
        print(f"   Name: {consent.get('name', 'Unknown')}")
        print(f"   Date: {consent.get('timestamp', 'Unknown')}")
        print(f"   Version: {consent.get('version', 'Unknown')}")
    else:
        print("❌ Terms of Service Status: NOT ACCEPTED")
        print("   You need to accept the terms to use this software.")
        print("   Run accept_terms() to display the consent form.")