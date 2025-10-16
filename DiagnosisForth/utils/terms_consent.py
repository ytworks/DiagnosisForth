import json
import socket
import threading
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, Any
import requests


def get_user_ip() -> str:
    """Get the user's IP address."""
    try:
        # Try to get external IP
        response = requests.get('https://api.ipify.org?format=json', timeout=5)
        if response.status_code == 200:
            return response.json()['ip']
    except:
        pass
    
    # Fallback to local IP
    try:
        hostname = socket.gethostname()
        local_ip = socket.gethostbyname(hostname)
        return local_ip
    except:
        return "Unknown"


def send_consent_to_slack(consent_data: Dict[str, Any]) -> None:
    """Send consent information to Slack."""
    try:
        base = "https://hooks.slack.com/services/"
        team = "T022ELD0Z33"
        channel = "B09LVFJPE2H"
        token = "jx7b0L0wvRqRnL6xbfzOTvpq"
        webhook_url = f"{base}{team}/{channel}/{token}"
        
        text = (
            f"📋 Terms of Service Accepted\n"
            f"━━━━━━━━━━━━━━━━━━━━━━\n"
            f"Name: {consent_data['name']}\n"
            f"IP Address: {consent_data['ip_address']}\n"
            f"Timestamp: {consent_data['timestamp']}\n"
            f"Version: {consent_data['version']}\n"
            f"Terms URL: https://github.com/ytworks/DiagnosisForth/blob/main/TERMS_OF_SERVICE.md\n"
            f"━━━━━━━━━━━━━━━━━━━━━━"
        )
        
        payload = {"text": text}
        thread = threading.Thread(
            target=lambda: requests.post(webhook_url, json=payload, timeout=5)
        )
        thread.daemon = True
        thread.start()
    except:
        pass


def check_existing_consent() -> Optional[Dict[str, Any]]:
    """Check if user has already consented to terms."""
    consent_file = Path.home() / ".diagnosisforth" / "consent.json"
    
    if consent_file.exists():
        try:
            with open(consent_file, 'r') as f:
                data = json.load(f)
                # Check if consent is for current version
                if data.get('version') == '1.0.0':
                    return data
        except:
            pass
    
    return None


def save_consent(consent_data: Dict[str, Any]) -> None:
    """Save consent data locally."""
    consent_dir = Path.home() / ".diagnosisforth"
    consent_dir.mkdir(exist_ok=True)
    consent_file = consent_dir / "consent.json"
    
    with open(consent_file, 'w') as f:
        json.dump(consent_data, f, indent=2)


def display_terms() -> None:
    """Display the terms of service."""
    print("\n" + "="*70)
    print(" TERMS OF SERVICE - IMPORTANT NOTICE ")
    print("="*70)
    print("""
This software includes functionality to monitor and log usage patterns.
By using this software, you consent to:

1. Collection of usage data for improvement and debugging
2. Monitoring of feature usage and performance metrics
3. Error tracking and reporting

DISCLAIMER OF LIABILITY:
━━━━━━━━━━━━━━━━━━━━━━━
THE AUTHOR(S) ARE COMPLETELY EXEMPT FROM ALL LIABILITY.
You use this software ENTIRELY AT YOUR OWN RISK.
The author(s) bear NO RESPONSIBILITY for any consequences,
damages, or issues arising from its use.

By accepting these terms, you:
• Acknowledge complete exemption of author liability
• Accept all risks associated with using this software
• Waive all rights to pursue claims against the author(s)
• Agree to indemnify the author(s) from any claims

Full terms available at: 
https://github.com/ytworks/DiagnosisForth/blob/main/TERMS_OF_SERVICE.md
""")
    print("="*70)


def request_consent() -> bool:
    """Request explicit consent from the user."""
    display_terms()
    
    print("\nTo use this software, you must explicitly accept the terms above.")
    print("Your acceptance will be logged for compliance purposes.\n")
    
    # Get user's name/signature
    name = input("Please enter your full name as digital signature: ").strip()
    if not name:
        print("❌ Name is required for consent.")
        return False
    
    # Confirm acceptance
    response = input(f"\n{name}, do you accept these terms? (yes/no): ").strip().lower()
    
    if response == 'yes':
        # Collect consent data
        consent_data = {
            'name': name,
            'ip_address': get_user_ip(),
            'timestamp': datetime.now().isoformat(),
            'version': '1.0.0',
            'accepted': True
        }
        
        # Save locally
        save_consent(consent_data)
        
        # Send to Slack (in background, non-blocking)
        send_consent_to_slack(consent_data)
        
        print("\n✅ Thank you for accepting the terms of service.")
        print("You may now use the software.\n")
        return True
    else:
        print("\n❌ You must accept the terms to use this software.")
        print("Exiting...\n")
        return False


def is_jupyter() -> bool:
    """Check if code is running in Jupyter environment."""
    try:
        from IPython import get_ipython
        if get_ipython() is not None:
            return True
    except ImportError:
        pass
    return False


def request_consent_jupyter():
    """Request consent in Jupyter notebook environment - automatically displays form."""
    try:
        from IPython.display import display, HTML, Javascript
        import ipywidgets as widgets
        import time
        
        # Display terms
        terms_html = """
        <div style="border: 2px solid #ff6b6b; padding: 20px; margin: 20px 0; background: #fff5f5;">
            <h2 style="color: #c92a2a;">⚠️ TERMS OF SERVICE - IMPORTANT NOTICE</h2>
            
            <p><strong>This software includes functionality to monitor and log usage patterns.</strong></p>
            
            <h3>By using this software, you consent to:</h3>
            <ul>
                <li>Collection of usage data for improvement and debugging</li>
                <li>Monitoring of feature usage and performance metrics</li>
                <li>Error tracking and reporting</li>
            </ul>
            
            <div style="background: #ffe0e0; padding: 15px; margin: 15px 0; border-left: 4px solid #ff6b6b;">
                <h3 style="color: #a61e1e;">⚠️ DISCLAIMER OF LIABILITY</h3>
                <p><strong>THE AUTHOR(S) ARE COMPLETELY EXEMPT FROM ALL LIABILITY.</strong></p>
                <p>You use this software ENTIRELY AT YOUR OWN RISK.</p>
                <p>The author(s) bear NO RESPONSIBILITY for any consequences, damages, or issues arising from its use.</p>
            </div>
            
            <h3>By accepting these terms, you:</h3>
            <ul>
                <li>Acknowledge complete exemption of author liability</li>
                <li>Accept all risks associated with using this software</li>
                <li>Waive all rights to pursue claims against the author(s)</li>
                <li>Agree to indemnify the author(s) from any claims</li>
            </ul>
            
            <p><em>Full terms available at: <a href="https://github.com/ytworks/DiagnosisForth/blob/main/TERMS_OF_SERVICE.md" target="_blank">GitHub - TERMS_OF_SERVICE.md</a></em></p>
        </div>
        """
        
        display(HTML(terms_html))
        
        # Create input widgets
        print("\n" + "="*70)
        print("To use this software, you must explicitly accept the terms above.")
        print("Your acceptance will be logged for compliance purposes.")
        print("="*70 + "\n")
        
        name_widget = widgets.Text(
            placeholder='Enter your full name',
            description='Digital Signature:',
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='500px')
        )
        
        accept_button = widgets.Button(
            description='I Accept the Terms',
            button_style='success',
            icon='check'
        )
        
        decline_button = widgets.Button(
            description='I Decline',
            button_style='danger',
            icon='times'
        )
        
        output = widgets.Output()
        
        # Flag to track if consent was given
        consent_given = {'value': False}
        
        def on_accept(b):
            with output:
                output.clear_output()
                name = name_widget.value.strip()
                
                if not name:
                    print("❌ Name is required for consent.")
                    return
                
                # Collect consent data
                consent_data = {
                    'name': name,
                    'ip_address': get_user_ip(),
                    'timestamp': datetime.now().isoformat(),
                    'version': '1.0.0',
                    'accepted': True,
                    'environment': 'jupyter'
                }
                
                # Save locally
                save_consent(consent_data)
                
                # Send to Slack
                send_consent_to_slack(consent_data)
                
                print("✅ Thank you for accepting the terms of service.")
                print("The module will be available after this cell completes.")
                
                # Hide widgets after acceptance
                name_widget.layout.display = 'none'
                accept_button.layout.display = 'none'
                decline_button.layout.display = 'none'
                
                consent_given['value'] = True
        
        def on_decline(b):
            with output:
                output.clear_output()
                print("❌ You must accept the terms to use this software.")
                print("The module will not be available until you accept the terms.")
        
        accept_button.on_click(on_accept)
        decline_button.on_click(on_decline)
        
        # Display widgets
        display(widgets.VBox([
            name_widget,
            widgets.HBox([accept_button, decline_button]),
            output
        ]))
        
        # Wait for user interaction (max 5 minutes)
        max_wait = 300  # 5 minutes
        start_time = time.time()
        
        while not consent_given['value'] and (time.time() - start_time) < max_wait:
            time.sleep(0.5)
        
        if consent_given['value']:
            return True
        else:
            return False
        
    except ImportError:
        print("⚠️ Jupyter widgets not available. Using text-based consent.")
        return request_consent()


def ensure_consent() -> bool:
    """
    Ensure user has consented to terms before using the software.
    Returns True if consent is given, False otherwise.
    """
    # Check for existing consent
    existing = check_existing_consent()
    
    if existing:
        # User has already consented for this version
        return True
    
    # Check if running in Jupyter
    if is_jupyter():
        return request_consent_jupyter()
    else:
        # Request new consent in terminal
        return request_consent()