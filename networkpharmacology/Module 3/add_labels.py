# /// script
# requires-python = ">=3.10"
# dependencies = ["pillow"]
# ///

from PIL import Image, ImageDraw, ImageFont
import math

def draw_arrow(draw, start, end, color="red", width=5, arrow_length=20, arrow_angle=math.pi/6):
    # Draw line
    draw.line([start, end], fill=color, width=width)
    
    # Calculate angle of the line
    angle = math.atan2(end[1] - start[1], end[0] - start[0])
    
    # Calculate coordinates for arrow head
    angle1 = angle + math.pi + arrow_angle
    angle2 = angle + math.pi - arrow_angle
    
    p1 = (end[0] + arrow_length * math.cos(angle1), end[1] + arrow_length * math.sin(angle1))
    p2 = (end[0] + arrow_length * math.cos(angle2), end[1] + arrow_length * math.sin(angle2))
    
    # Draw arrow head
    draw.polygon([end, p1, p2], fill=color)

# Load the image
img = Image.open("/home/hmalik342/Working/Python-BioPython/Demo Project/Module 3/CTSS_Pocket8.png")
draw = ImageDraw.Draw(img)

# Center of image (pocket)
center_x = img.width // 2
center_y = img.height // 2

# Text and arrow positions
text_pos = (150, 150)
arrow_start = (300, 200)
# Arrow end slightly offset from the exact center so it points TO the pocket, not over it
arrow_end = (center_x - 30, center_y - 30)

# Try to use a default truetype font, fallback to default bitmap font if not found
try:
    font = ImageFont.truetype("DejaVuSans-Bold.ttf", 40)
except IOError:
    font = ImageFont.load_default()

# Draw text
text = "Druggability Score: 0.713\n(Best Pocket)"
# simple shadow for text
draw.text((text_pos[0]+2, text_pos[1]+2), text, font=font, fill="gray")
draw.text(text_pos, text, font=font, fill="black")

# Draw arrow pointing to the pocket
draw_arrow(draw, arrow_start, arrow_end, color="#D32F2F", width=8, arrow_length=30)

# Save the final labeled image
img.save("/home/hmalik342/.gemini/antigravity/brain/eb3deb5f-076a-4c2f-aa4a-39ba3166e0e7/CTSS_Pocket8.png")
print("Image labeled and saved successfully.")
