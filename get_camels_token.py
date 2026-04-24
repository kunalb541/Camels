"""
get_camels_token.py — authenticate for CAMELS HTTPS download.

Run this once in your terminal:
    python get_camels_token.py

It opens a browser tab (or prints a URL), you authorize, paste the code,
and it saves an HTTPS access token to ~/.globus/camels_https_token.txt
for use by download_camels.py.
"""
import pathlib
import globus_sdk

CLI_CLIENT_ID   = "95fdeba8-fac2-42bd-a357-e068d82ff78e"
COLLECTION_ID   = "58bdcd24-6590-11ec-9b60-f9dfb1abb183"
HTTPS_SCOPE     = f"https://auth.globus.org/scopes/{COLLECTION_ID}/https"
TOKEN_FILE      = pathlib.Path.home() / ".globus" / "camels_https_token.txt"

client = globus_sdk.NativeAppAuthClient(CLI_CLIENT_ID)
client.oauth2_start_flow(requested_scopes=[HTTPS_SCOPE])

auth_url = client.oauth2_get_authorize_url()
print("\nOpen this URL in your browser to authorize:\n")
print(auth_url)
print()

auth_code = input("Paste the authorization code here: ").strip()

tokens = client.oauth2_exchange_code_for_tokens(auth_code)
https_token = tokens.by_resource_server[COLLECTION_ID]["access_token"]

TOKEN_FILE.parent.mkdir(parents=True, exist_ok=True)
TOKEN_FILE.write_text(https_token)
print(f"\nToken saved to {TOKEN_FILE}")
print("Now run:  python download_camels.py")
